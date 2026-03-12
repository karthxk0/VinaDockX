import os
import sys
import argparse
import glob
import platform
import subprocess
import shutil
from pathlib import Path

# Try importing required libraries
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("\n[CRITICAL ERROR] RDKit is not installed.")
    print("Please install it using: pip install rdkit")
    sys.exit(1)

try:
    from tqdm import tqdm
except ImportError:
    print("\n[CRITICAL ERROR] tqdm is not installed.")
    print("Please install it using: pip install tqdm")
    sys.exit(1)

# --------------------------------------------------------------------------------------
# BANNER
# --------------------------------------------------------------------------------------
def print_banner():
    banner = r"""
======================================================================                                                              
  _____                _      _       
 |  __ \              | |    (_)      
 | |__) | __ ___ _ __ | |     _  __ _ 
 |  ___/ '__/ _ \ '_ \| |    | |/ _` |
 | |   | | |  __/ |_) | |____| | (_| |
 |_|   |_|  \___| .__/|______|_|\__, |.py
                | |              __/ |
                |_|             |___/ 
                                                                         
  > Version 2.3  |  > Built by karthxk (https://karthxk0.github.io/)          
   
====================================================================== 
"""
    print(banner)

# --------------------------------------------------------------------------------------
# HELPER: Check if 3D
# --------------------------------------------------------------------------------------
def is_3d(sdf_file):
    """
    Check whether the first molecule in the SDF file is 3D or 2D.
    Returns True if any atom has a z-coordinate that is not close to zero.
    """
    try:
        with open(sdf_file, "r") as f:
            lines = f.readlines()
    except Exception as e:
        return False

    if len(lines) < 4:
        return False

    counts_line = lines[3].strip()
    tokens = counts_line.split()
    if not tokens:
        return False
    try:
        num_atoms = int(tokens[0])
    except ValueError:
        return False

    for i in range(4, 4 + num_atoms):
        if i >= len(lines):
            break
        parts = lines[i].split()
        if len(parts) < 3:
            continue
        try:
            z = float(parts[2])
        except ValueError:
            continue
        if abs(z) > 1e-3:
            return True
    return False

# --------------------------------------------------------------------------------------
# HELPER: Get Meeko Command
# --------------------------------------------------------------------------------------
def get_meeko_command():
    system = platform.system()
    if system == "Windows":
        cmd = "mk_prepare_ligand"
    else:
        cmd = "mk_prepare_ligand.py"
    
    from shutil import which
    if which(cmd.replace(".py", "")) is None and which(cmd) is None:
        return cmd
    return cmd

# --------------------------------------------------------------------------------------
# PROCESSING PIPELINE
# --------------------------------------------------------------------------------------
def process_ligands(input_files, output_dir, save_intermediates):
    
    meeko_cmd = get_meeko_command()
    
    # Stats Tracking
    stats = {
        "total_inputs": 0,
        "dist": {"sdf": 0, "mol2": 0},
        "initial_2d": 0,
        "initial_3d": 0,
        "converted_2d_to_3d": 0,
        "pdbqt_generated": 0,
        "success": 0,
        "failures": [] 
    }

    # Create Main Output Directory
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except Exception as e:
            print(f"[CRITICAL] Could not create output dir: {e}")
            sys.exit(1)

    # Paths for subfolders (DO NOT create yet)
    sdf_out_dir = os.path.join(output_dir, "Converted to SDF")
    threed_out_dir = os.path.join(output_dir, "Converted to 3D")
    
    print("\nStarting Processing Pipeline...")
    print("-" * 60)

    stats["total_inputs"] = len(input_files)
    
    # Initialize Progress Bar
    pbar = tqdm(input_files, unit="mol", desc="Processing")
    
    for file_path in pbar:
        filename = os.path.basename(file_path)
        name_no_ext = os.path.splitext(filename)[0]
        ext = os.path.splitext(filename)[1].lower()
        
        # Update Distribution
        if "sdf" in ext: stats["dist"]["sdf"] += 1
        elif "mol2" in ext: stats["dist"]["mol2"] += 1
        
        pbar.set_postfix_str(f"Current: {filename}")
        
        temp_files_to_delete = []
        
        try:
            # =================================================================
            # STAGE 1: FORMAT STANDARDIZATION (Convert MOL2 -> SDF)
            # =================================================================
            current_sdf_path = file_path
            
            if ext == '.mol2':
                tqdm.write(f"\n[>] {filename}: Converting to SDF...")
                
                # Try standard load
                mol = Chem.MolFromMol2File(file_path, removeHs=False)
                
                # Fallback: PyMOL-generated Mol2s often fail RDKit strict parsing
                if mol is None:
                    tqdm.write("    [!] Standard load failed, trying sanitize=False...")
                    mol = Chem.MolFromMol2File(file_path, removeHs=False, sanitize=False)
                    if mol:
                        mol.UpdatePropertyCache(strict=False)
                        try:
                            Chem.SanitizeMol(mol)
                        except:
                            tqdm.write("    [!] Sanitization warning (structure might be messy).")

                if mol is None:
                    raise ValueError(f"RDKit failed to load .MOL2 file. (Check format)")
                
                if save_intermediates:
                    if not os.path.exists(sdf_out_dir): os.makedirs(sdf_out_dir)
                    current_sdf_path = os.path.join(sdf_out_dir, f"{name_no_ext}.sdf")
                else:
                    current_sdf_path = os.path.join(output_dir, f"temp_{name_no_ext}_from_{ext[1:]}.sdf")
                    temp_files_to_delete.append(current_sdf_path)
                
                writer = Chem.SDWriter(current_sdf_path)
                writer.write(mol)
                writer.close()
            
            # =================================================================
            # STAGE 2: 3D CHECK & CONVERSION
            # =================================================================
            is_structure_3d = is_3d(current_sdf_path)
            path_for_meeko = current_sdf_path
            
            if not is_structure_3d:
                stats["initial_2d"] += 1
                tqdm.write(f"    - Detected 2D geometry. Generating 3D coordinates...")
                
                mol = Chem.SDMolSupplier(current_sdf_path, removeHs=False)[0]
                if mol is None: raise ValueError("Could not read SDF for 3D generation.")
                
                mol = Chem.AddHs(mol)
                res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                if res == -1:
                    res = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                    if res == -1: raise ValueError("3D Embedding failed.")
                
                AllChem.MMFFOptimizeMolecule(mol)
                
                if save_intermediates:
                    if not os.path.exists(threed_out_dir): os.makedirs(threed_out_dir)
                    path_for_meeko = os.path.join(threed_out_dir, f"{name_no_ext}_3d.sdf")
                else:
                    path_for_meeko = os.path.join(output_dir, f"temp_{name_no_ext}_3d.sdf")
                    temp_files_to_delete.append(path_for_meeko)
                
                writer = Chem.SDWriter(path_for_meeko)
                writer.write(mol)
                writer.close()
                stats["converted_2d_to_3d"] += 1
                
            else:
                stats["initial_3d"] += 1
                tqdm.write(f"    - Detected 3D geometry. Refining (adding Hs/Minimizing)...")
                
                mol = Chem.SDMolSupplier(current_sdf_path, removeHs=False)[0]
                mol = Chem.AddHs(mol, addCoords=True)
                AllChem.MMFFOptimizeMolecule(mol)
                
                # Always write to a new file for Meeko to ensure Hs are saved
                path_for_meeko = os.path.join(output_dir, f"temp_{name_no_ext}_refined.sdf")
                temp_files_to_delete.append(path_for_meeko)
                
                writer = Chem.SDWriter(path_for_meeko)
                writer.write(mol)
                writer.close()

            # =================================================================
            # STAGE 3: MEEKO PDBQT GENERATION
            # =================================================================
            pdbqt_name = f"{name_no_ext}_PrepLig.pdbqt"
            pdbqt_path = os.path.join(output_dir, pdbqt_name)
            
            # --- CRITICAL FIX FOR WINDOWS PATHS WITH SPACES ---
            # We explicitly quote paths and pass a string to shell on Windows
            if platform.system() == "Windows":
                cmd_string = f'{meeko_cmd} -i "{path_for_meeko}" -o "{pdbqt_path}"'
                result = subprocess.run(cmd_string, capture_output=True, text=True, shell=True)
            else:
                # Linux/Mac usually handles lists correctly
                cmd_args = [meeko_cmd, "-i", path_for_meeko, "-o", pdbqt_path]
                result = subprocess.run(cmd_args, capture_output=True, text=True, shell=False)
            
            if result.returncode == 0 and os.path.exists(pdbqt_path):
                stats["pdbqt_generated"] += 1
                stats["success"] += 1
            else:
                # Clean error message
                err_msg = result.stderr.strip() if result.stderr else "Unknown Meeko Error"
                raise RuntimeError(f"Meeko failed. Stderr: {err_msg}")

        except Exception as e:
            stats["failures"].append((filename, str(e)))
        
        finally:
            # Clean up temp files
            for temp_f in temp_files_to_delete:
                if os.path.exists(temp_f):
                    try: os.remove(temp_f)
                    except: pass

    pbar.close()
    
    # ----------------------------------------------------------------------------------
    # FINAL DETAILED SUMMARY
    # ----------------------------------------------------------------------------------
    print("\n" + "="*80)
    print(f"{'PROCESSING SUMMARY':^80}")
    print("="*80)
    
    print(f"\n1. INPUT FILES")
    print(f"   Total Input Files    : {stats['total_inputs']}")
    print(f"   Distribution         : SDF: {stats['dist']['sdf']} | MOL2: {stats['dist']['mol2']}")
    
    print(f"\n2. GEOMETRY ANALYSIS (Initial State)")
    print(f"   2D Files Detected    : {stats['initial_2d']}")
    print(f"   3D Files Detected    : {stats['initial_3d']}")
    
    print(f"\n3. CONVERSIONS & GENERATION")
    print(f"   2D -> 3D Converted   : {stats['converted_2d_to_3d']} (Saved in 'Converted to 3D' if enabled)")
    print(f"   PDBQT Files Created  : {stats['pdbqt_generated']}")
    
    print(f"\n4. OVERALL STATUS")
    print(f"   Successful Molecules : {stats['success']}")
    print(f"   Failed Molecules     : {len(stats['failures'])}")
    
    if stats['failures']:
        print("\n" + "-"*80)
        print("ERROR LOG")
        print("-"*80)
        print(f"{'Filename':<30} | {'Error Message'}")
        print("-"*80)
        for fname, err in stats['failures']:
            print(f"{fname:<30} | {err}")
    
    print("\n" + "="*80 + "\n")

# --------------------------------------------------------------------------------------
# INPUT PARSING
# --------------------------------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(description="Prepare ligands for AutoDock Vina.")
    parser.add_argument("-i", "--input", help="Input files (comma separated) or directories")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-s", "--save_intermediates", type=int, choices=[0, 1], 
                        help="Save intermediate files? 1=Yes, 0=No")
    return parser.parse_args()

def get_files_from_input(input_str):
    allowed_ext = ('.sdf', '.mol2')
    files_found = []
    
    raw_inputs = [x.strip().strip('"').strip("'") for x in input_str.split(',')]
    
    for item in raw_inputs:
        path_item = Path(item)
        if path_item.is_dir():
            for ext in allowed_ext:
                found = list(path_item.rglob(f"*{ext}"))
                files_found.extend([str(p) for p in found])
        elif path_item.is_file():
            if str(path_item).lower().endswith(allowed_ext):
                files_found.append(str(path_item))
        else:
            glob_found = glob.glob(item)
            for g in glob_found:
                if g.lower().endswith(allowed_ext):
                    files_found.append(g)
                    
    return list(set(files_found))

# --------------------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------------------
def main():
    print_banner()
    
    args = parse_arguments()
    
    if len(sys.argv) == 1:
        # Interactive Mode
        print("INTERACTIVE MODE\n")
        
        input_str = input(">> Please enter Ligand files (SDF/MOL2) or Folder paths (separated by comma): ")
        input_files = get_files_from_input(input_str)
        
        if not input_files:
            print("[!] No valid files found. Exiting.")
            return

        print(f"   Found {len(input_files)} valid files.")

        try:
            save_int_input = input(">> Save generated SDF/3D files separately? (1: Yes, 0: No) : ")
            if save_int_input.strip() == "":
                save_intermediates = False
            else:
                save_intermediates = bool(int(save_int_input))
        except ValueError:
            save_intermediates = False
            print("   Invalid input. Defaulting to No.")

        output_dir = input(">> Enter Output Directory: ").strip().strip('"').strip("'")
        if not output_dir:
            print("[!] Output directory cannot be empty.")
            return

    else:
        # CLI Mode
        if not args.input or not args.output:
            print("[!] CLI Usage Error: -i and -o are required in CLI mode.")
            print("    Usage: python PrepLig.py -i <files> -o <dir> [-s 0/1]")
            return
        
        input_files = get_files_from_input(args.input)
        if not input_files:
            print("[!] No valid files found based on input arguments.")
            return
            
        save_intermediates = bool(args.save_intermediates) if args.save_intermediates is not None else False
        output_dir = args.output

    process_ligands(input_files, output_dir, save_intermediates)

if __name__ == "__main__":
    main()