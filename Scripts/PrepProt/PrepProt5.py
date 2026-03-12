#!/usr/bin/env python3
import sys
import os
import argparse
import subprocess
import platform
import importlib.util
from tqdm import tqdm

# --- Configuration ---
VERSION = "2.5"

# --- Helper Functions ---

def print_banner():
    banner = r"""
===========================================================
  _____                _____           _    
 |  __ \              |  __ \         | |   
 | |__) | __ ___ _ __ | |__) | __ ___ | |_  
 |  ___/ '__/ _ \ '_ \|  ___/ '__/ _ \| __| 
 | |   | | |  __/ |_) | |   | | | (_) | |_  
 |_|   |_|  \___| .__/|_|   |_|  \___/ \__|.py 
                | |                         
                |_|                         
                                                                                                                                                              
    > Version: 2.5                                  
    > Designed by: karthxk (https://karthxk0.github.io/)   
      
===========================================================
    """
    print(banner)

def check_dependencies():
    """Silently checks dependencies, only alerting on failure."""
    missing = []
    packages = ['meeko', 'tqdm', 'Bio'] 
    for pkg in packages:
        if importlib.util.find_spec(pkg) is None:
            missing.append(pkg)
    
    if missing:
        print(f"\n[!] Missing Dependencies: {', '.join(missing)}")
        print("    Run: pip install meeko tqdm biopython")
        sys.exit(1)

def normalize_path(path_str):
    if not path_str: return ""
    return path_str.strip().strip('"').strip("'")

def parse_flexible_residues(flex_input):
    """
    Parses flexible residues input.
    Accepts: A:159, B:241 OR [('A', 159), ('B', 241)]
    Returns: A:159,B:241
    """
    if not flex_input: return None
    s = flex_input.strip()
    
    # Handle Python List format
    if '[' in s and ']' in s:
        import re
        matches = re.findall(r"['\"]([a-zA-Z0-9]+)['\"],\s*(\d+)", s)
        if matches:
            return ",".join([f"{c}:{r}" for c, r in matches])

    # Handle Standard format
    return s.replace(" ", "")

def get_all_pdb_files(input_str):
    inputs = [normalize_path(i) for i in input_str.split(',')]
    pdb_files = []
    for item in inputs:
        if os.path.isfile(item) and item.lower().endswith('.pdb'):
            pdb_files.append(os.path.abspath(item))
        elif os.path.isdir(item):
            for root, _, files in os.walk(item):
                for file in files:
                    if file.lower().endswith('.pdb'):
                        pdb_files.append(os.path.abspath(os.path.join(root, file)))
    return sorted(list(set(pdb_files)))

def clean_pdb_biopython(input_pdb, output_pdb):
    """
    Uses BioPython to strip everything except Protein Heavy Atoms.
    Removes Water, Ligands, Ions, and ALL Hydrogens (Polar & Non-Polar).
    """
    try:
        from Bio import PDB
        
        class CleanSelect(PDB.Select):
            def accept_residue(self, residue):
                # Keep standard amino acids only
                return PDB.is_aa(residue, standard=True)

            def accept_atom(self, atom):
                # Remove ALL hydrogens to ensure clean slate for Meeko
                if atom.element.upper() == 'H': return False
                # Double check for water/ions
                if atom.parent.resname in ['HOH', 'WAT']: return False
                return True

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", input_pdb)
        
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_pdb, select=CleanSelect())
        return True
        
    except Exception as e:
        return False

def run_meeko(input_pdb, output_pdbqt, flex_residues=None):
    """Runs mk_prepare_receptor."""
    system_os = platform.system()
    cmd_base = "mk_prepare_receptor" if system_os == "Windows" else "mk_prepare_receptor.py"

    cmd = [cmd_base, "--read_pdb", input_pdb, "-o", output_pdbqt, "-p", "-a"]
    if flex_residues:
        cmd.extend(["-f", flex_residues])

    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True, "Success"
    except subprocess.CalledProcessError as e:
        err = e.stderr
        if "Explicit valence" in err:
            return False, "Geometry Error (Valence)"
        if "ValueError" in err:
             return False, "Meeko Value Error"
        return False, err.strip()
    except FileNotFoundError:
        return False, f"Missing Exe: {cmd_base}"

# --- Main Interface ---

def main():
    check_dependencies()
    print_banner()

    parser = argparse.ArgumentParser(description="PrepProt.py")
    parser.add_argument("-i", "--input", help="Input PDBs/Folders")
    parser.add_argument("-f", "--flex", help="Flexible Residues")
    parser.add_argument("-o", "--outdir", help="Output Directory")
    parser.add_argument("-s", "--subfolders", action="store_true", help="Use Subfolders")
    args = parser.parse_args()

    # --- Input Collection ---
    if len(sys.argv) > 1:
        # CLI Mode
        input_raw = args.input
        flex_raw = args.flex
        out_dir = normalize_path(args.outdir)
        use_subfolders = args.subfolders
        print(f" [Mode] CLI Execution")
    else:
        # Interactive Mode (Aesthetic)
        print(" [Configuration Step]")
        input_raw = input("   Input PDB Files/Folders (separated by comma): ")
        
        # Compact Logic for Yes/No
        ask_flex = input("   Define Flexible Residues? [1=Yes, 0=No]: ") or "0"
        flex_raw = None
        if ask_flex.strip() == "1":
            print("      (Sample format: A:159, B:241 or [('A', 159), ('B', 241)])")
            flex_raw = input("      > Residues: ")
        
        out_dir = normalize_path(input("  [?] Output Directory: "))
        
        ask_sub = input("   Organize each PDB into Sub-directory? (1=Yes, 0=No): ") or "0"
        use_subfolders = (ask_sub.strip() == "1")
        print("") # Spacer

    # --- Processing ---
    if not input_raw:
        print("\n [!] Error: No input provided."); sys.exit(1)

    pdb_files = get_all_pdb_files(input_raw)
    flex_str = parse_flexible_residues(flex_raw)
    
    if not pdb_files:
        print("\n [!] Error: No PDB files found."); sys.exit(1)
    
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    print(f" [Processing Queue] {len(pdb_files)} Proteins")
    print("-" * 60)

    failed_log = []

    # Clean Progress Bar
    # 'ncols=80' keeps it from wrapping on small screens
    # 'bar_format' simplifies the visual
    pbar_fmt = "{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]"
    
    with tqdm(total=len(pdb_files), bar_format=pbar_fmt, ncols=75) as pbar:
        for pdb_path in pdb_files:
            base = os.path.splitext(os.path.basename(pdb_path))[0]
            pbar.set_description(f" > {base[:10]:<10}")
            
            target_dir = os.path.join(out_dir, base) if use_subfolders else out_dir
            if use_subfolders: os.makedirs(target_dir, exist_ok=True)
            
            temp_pdb = os.path.join(target_dir, f"{base}_temp.pdb")
            final_qt = os.path.join(target_dir, f"{base}_PrepProt")
            
            # 1. Clean (BioPython)
            if clean_pdb_biopython(pdb_path, temp_pdb):
                # 2. Convert (Meeko)
                ok, msg = run_meeko(temp_pdb, final_qt, flex_str)
                if not ok:
                    failed_log.append((base, msg))
                if os.path.exists(temp_pdb): os.remove(temp_pdb)
            else:
                failed_log.append((base, "Cleaning Failed"))
            
            pbar.update(1)

    # --- Final Summary ---
    success_count = len(pdb_files) - len(failed_log)
    
    print("\n" + "="*60)
    print(f"{'Summary Report':^60}")
    print("-" * 60)
    print(f"  Total Processed : {len(pdb_files)}")
    print(f"  Successful      : {success_count}")
    print(f"  Failed          : {len(failed_log)}")
    
    if failed_log:
        print("-" * 60)
        print("  [!] Failed Molecules:")
        for name, reason in failed_log:
            # Truncate long errors for aesthetics
            clean_reason = reason.replace("\n", " ")
            if len(clean_reason) > 40: clean_reason = clean_reason[:37] + "..."
            print(f"   - {name:<15} : {clean_reason}")
            
    print("="*60 + "\n")

if __name__ == "__main__":
    main()