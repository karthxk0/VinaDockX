import os
import sys
import argparse
import re

# ==========================================
#              DEPENDENCY CHECK
# ==========================================
def check_dependencies():
    """Checks if required packages are installed."""
    missing = []
    try:
        import numpy as np
    except ImportError:
        missing.append("numpy")
    
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        missing.append("biopython")
        
    if missing:
        print("\n[!] CRITICAL ERROR: Missing required packages.")
        print(f"    Please install: {', '.join(missing)}")
        print("    Run: pip install " + " ".join(missing))
        sys.exit(1)

# Perform check before importing logic dependent on them
check_dependencies()
import numpy as np
from Bio.PDB import PDBParser

# ==========================================
#              HELPER FUNCTIONS
# ==========================================

def print_banner():
    """Prints the aesthetic banner."""
    banner = r"""
========================================================================
    
   _____      _     _  _____             __ _        _____            
  / ____|    (_)   | |/ ____|           / _(_)      / ____|           
 | |  __ _ __ _  __| | |     ___  _ __ | |_ _  __ _| |  __  ___ _ __  
 | | |_ | '__| |/ _` | |    / _ \| '_ \|  _| |/ _` | | |_ |/ _ \ '_ \ 
 | |__| | |  | | (_| | |___| (_) | | | | | | | (_| | |__| |  __/ | | |
  \_____|_|  |_|\__,_|\_____\___/|_| |_|_| |_|\__, |\_____|\___|_| |_|.py
                                               __/ |                  
                                              |___/                   
                                              
    > Version: 2.3                                     
    > Designed by: karthxk (https://karthxk0.github.io/)   
       
========================================================================
    """
    print(banner)

def clean_path(path_str):
    """
    Removes surrounding quotes from file paths and normalizes slashes.
    """
    if not path_str:
        return ""
    clean = path_str.strip().strip('"').strip("'")
    return os.path.normpath(clean)

def extract_residues(text):
    """
    Uses Regex to find residues in two formats:
    1. A:12, B:15
    2. [('A', 12), ('B', 15)]
    """
    residues = []
    
    # Pattern 1: A:12 (Matches Letter, optional space, colon, optional space, digits)
    pattern1 = r'([A-Za-z])\s*:\s*(\d+)'
    for chain, resnum in re.findall(pattern1, text):
        residues.append((chain.upper(), int(resnum)))
        
    # Pattern 2: ('A', 12) or ("A", 12)
    pattern2 = r"\(\s*['\"]([A-Za-z])['\"]\s*,\s*(\d+)\s*\)"
    for chain, resnum in re.findall(pattern2, text):
        residues.append((chain.upper(), int(resnum)))
        
    # Remove duplicates while preserving order
    unique_res = []
    for res in residues:
        if res not in unique_res:
            unique_res.append(res)
            
    return unique_res

def parse_residues_input(input_val):
    """
    Parses residue input by reading the file or string and extracting matching patterns.
    """
    input_val = clean_path(input_val)
    content = ""
    
    # Check if input is a file
    if os.path.isfile(input_val):
        try:
            with open(input_val, 'r') as f:
                content = f.read()
        except Exception as e:
            raise ValueError(f"Error reading residue file: {e}")
    else:
        # Treat as raw string input
        content = input_val

    residues = extract_residues(content)
    
    if not residues:
        raise ValueError("No valid residues found. Ensure format is A:12 or [('A', 12)].")
        
    return residues

def format_vector_string(vec):
    """Formats a numpy vector to a string like [x, y, z] with commas."""
    return f"[{vec[0]:.3f}, {vec[1]:.3f}, {vec[2]:.3f}]"

# ==========================================
#              CORE LOGIC
# ==========================================

def calculate_protein_dimensions(pdb_file, padding=5.0):
    """Calculates center and size for blind docking (whole protein)."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_file)
    except Exception as e:
        raise ValueError(f"Failed to parse PDB file: {e}")

    atom_coords = np.array([atom.coord for atom in structure.get_atoms()])
    
    if len(atom_coords) == 0:
        raise ValueError("No atoms found in the PDB structure.")

    min_coord = np.min(atom_coords, axis=0)
    max_coord = np.max(atom_coords, axis=0)
    
    center = (min_coord + max_coord) / 2
    size = (max_coord - min_coord) + padding
    return center, size

def calculate_binding_site_center(pdb_file, residues, padding=5.0):
    """Calculates center and size for targeted docking based on residues."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_file)
    except Exception as e:
        raise ValueError(f"Failed to parse PDB file: {e}")

    model = structure[0]
    coords = []
    
    for chain_id, res_num in residues:
        try:
            res_num = int(res_num)
            chain = model[chain_id]
            residue = chain[res_num]
            if 'CA' in residue:
                coords.append(residue['CA'].coord)
            else:
                coords.append(list(residue.get_atoms())[0].coord)
        except KeyError:
            print(f"Warning: Residue {res_num} in chain {chain_id} not found.")
            continue
        except Exception as e:
            print(f"Error processing residue {chain_id}:{res_num} - {e}")
            continue

    if len(coords) == 0:
        raise ValueError("No valid residues found to calculate grid.")
    
    center = np.mean(coords, axis=0)
    size = np.ptp(coords, axis=0) + padding
    return center, size

def generate_file_content(center, size, pdb_source_path, mode_type, is_config=False):
    """Generates the string content for the output file."""
    
    if mode_type == 'blind':
        mode_desc = "Blind (whole protein)"
    else:
        mode_desc = "Targeted (based on residues provided)"

    content = (
        f"# Generated by GridnConfigGen.py v2.3\n"
        f"# Designed by karthxk (https://karthxk0.github.io/)\n\n"
        f"# Input: {pdb_source_path}\n"
        f"# Mode: {mode_desc}\n\n"
        f"# Grid Box Center\n"
        f"center_x = {center[0]:.3f}\n"
        f"center_y = {center[1]:.3f}\n"
        f"center_z = {center[2]:.3f}\n\n"
        f"# Grid Box Size\n"
        f"size_x = {size[0]:.3f}\n"
        f"size_y = {size[1]:.3f}\n"
        f"size_z = {size[2]:.3f}\n"
    )

    if is_config:
        content += (
            f"\n# Docking Parameters\n"
            f"exhaustiveness = 32\n"
            f"verbosity = 1\n"
        )
        
    return content

def save_files(center, size, output_dir, pdb_source_path, mode_type, output_choice):
    """
    Handles file generation based on output_choice.
    1: Grid only
    2: Config only
    3: Both
    """
    
    pdb_basename = os.path.basename(pdb_source_path)
    pdb_name_no_ext = os.path.splitext(pdb_basename)[0]
    suffix = "Blind" if mode_type == 'blind' else "Targeted"

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created directory: {output_dir}")
        except Exception as e:
            print(f"Error creating directory {output_dir}: {e}")
            return

    generated_files = []

    # --- Generate Grid File (Option 1 or 3) ---
    if output_choice in [1, 3]:
        filename_grid = f"{pdb_name_no_ext}_Grid_{suffix}.txt"
        path_grid = os.path.join(output_dir, filename_grid)
        content_grid = generate_file_content(center, size, pdb_source_path, mode_type, is_config=False)
        
        try:
            with open(path_grid, "w") as f:
                f.write(content_grid)
            generated_files.append(filename_grid)
        except Exception as e:
            print(f"Error writing grid file: {e}")

    # --- Generate Config File (Option 2 or 3) ---
    if output_choice in [2, 3]:
        filename_conf = f"{pdb_name_no_ext}_Config_{suffix}.txt"
        path_conf = os.path.join(output_dir, filename_conf)
        content_conf = generate_file_content(center, size, pdb_source_path, mode_type, is_config=True)
        
        try:
            with open(path_conf, "w") as f:
                f.write(content_conf)
            generated_files.append(filename_conf)
        except Exception as e:
            print(f"Error writing config file: {e}")

    # --- Print Summary Message ---
    if len(generated_files) == 1:
        print(f"Generated Grid file {generated_files[0]} saved to {output_dir}")
    elif len(generated_files) == 2:
        print(f"Generated Grid file ({generated_files[0]}) and Config file ({generated_files[1]}) saved to {output_dir}")

# ==========================================
#              INTERACTION MODES
# ==========================================

def run_interactive():
    """Interactive Wizard Mode"""
    # 1. Select Mode
    print("Select the type of grid you would like to generate:")
    print("1 - Blind Docking (whole protein)")
    print("2 - Targeted Docking (specific binding site)")
    choice = input("Enter 1 or 2: ").strip()

    # 2. Select PDB
    pdb_file = input("Enter path to PDB file: ")
    pdb_file = clean_path(pdb_file)

    if not os.path.exists(pdb_file):
        print("Error: PDB file does not exist.")
        return

    # Check for residues input immediately if targeted mode is selected
    residues = None
    if choice == '2':
        res_input = input("Enter residues list OR path to .txt file (e.g. [('A', 12), ('B', 15)] or A:12, B:15 ): ")
        try:
            residues = parse_residues_input(res_input)
            print(f"Detected residues: {residues}")
        except ValueError as e:
            print(e)
            return

    # 3. Select Output Format
    print("\nSelect output format:")
    print("1. Grid file")
    print("2. Config file (for VinaDock)")
    print("3. Both Grid & Config")
    print("4. None (Just display grid parameters)")
    
    out_choice_str = input("Enter number (1-4): ").strip()
    if out_choice_str not in ['1', '2', '3', '4']:
        print("Invalid output choice.")
        return
    out_choice = int(out_choice_str)

    # 4. Select Output Directory (Skipped if Option 4)
    output_dir = ""
    if out_choice != 4:
        output_dir = input("\nEnter directory to save output file: ")
        output_dir = clean_path(output_dir)

    # 5. Perform Calculations
    center, size = None, None
    mode_type = ""

    # Hardcoding padding for interactive mode per user request
    padding_val = 5.0 

    try:
        if choice == '1':
            mode_type = "blind"
            center, size = calculate_protein_dimensions(pdb_file, padding=padding_val)
        elif choice == '2':
            mode_type = "targeted"
            center, size = calculate_binding_site_center(pdb_file, residues, padding=padding_val)
        else:
            print("Invalid choice.")
            return
    except Exception as e:
        print(f"Calculation Error: {e}")
        return

    # 6. Display Grid
    print(f"\nCalculated Grid Center: {format_vector_string(center)}")
    print(f"Calculated Grid Size:   {format_vector_string(size)}\n")

    # 7. Generate Files
    if out_choice == 4:
        print("No files generated.")
    else:
        save_files(center, size, output_dir, pdb_file, mode_type, out_choice)

def run_cli():
    """Command Line Interface Mode for Pipelines"""
    parser = argparse.ArgumentParser(description="GridnConfigGen v2.3 CLI")
    
    parser.add_argument("--mode", choices=['blind', 'targeted'], required=True, help="Docking mode")
    parser.add_argument("--pdb", required=True, help="Path to input PDB file")
    parser.add_argument("--type", choices=['grid', 'config', 'both', 'none'], default='grid', help="Output file type (default: grid)")
    parser.add_argument("--out", help="Directory to save output (Required unless --type is none)")
    parser.add_argument("--residues", help="List string or path to .txt file (Required for targeted mode)")
    parser.add_argument("--padding", type=float, default=5.0, help="Padding for bounding box (default 5.0)")

    args = parser.parse_args()
    
    if args.type != 'none' and not args.out:
        parser.error("--out is required unless --type is 'none'")

    pdb_file = clean_path(args.pdb)
    output_dir = clean_path(args.out) if args.out else ""

    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found.")
        sys.exit(1)

    try:
        if args.mode == 'blind':
            center, size = calculate_protein_dimensions(pdb_file, padding=args.padding)
        elif args.mode == 'targeted':
            if not args.residues:
                print("Error: --residues argument is required for targeted docking.")
                sys.exit(1)
            residues = parse_residues_input(args.residues)
            print(f"Detected residues: {residues}")
            center, size = calculate_binding_site_center(pdb_file, residues, padding=args.padding)
        
        print(f"Processing {os.path.basename(pdb_file)} in {args.mode} mode...")
        print(f"Calculated Grid Center: {format_vector_string(center)}")
        print(f"Calculated Grid Size:   {format_vector_string(size)}")

        if args.type == 'none':
            print("No files generated.")
        else:
            choice_map = {'grid': 1, 'config': 2, 'both': 3}
            save_files(center, size, output_dir, pdb_file, args.mode, choice_map[args.type])
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    print_banner()
    if len(sys.argv) > 1:
        run_cli()
    else:
        run_interactive()

if __name__ == "__main__":
    main()