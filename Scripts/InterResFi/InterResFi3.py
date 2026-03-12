import MDAnalysis as mda
import warnings
import os
import sys
import argparse
from pathlib import Path

# Suppress MDAnalysis warnings
warnings.filterwarnings('ignore')

# --- CONFIGURATION ---
SCRIPT_NAME = "FindInterRes"
VERSION = "2.3"
AUTHOR = "karthxk"
URL = "https://karthxk.github.io/"

def print_header():
    """Prints the stylized ASCII header."""
    print(r"""
 ==============================================================
  ______ _           _ _____       _            _____           
 |  ____(_)         | |_   _|     | |          |  __ \          
 | |__   _ _ __   __| | | |  _ __ | |_ ___ _ __| |__) |___  ___ 
 |  __| | | '_ \ / _` | | | | '_ \| __/ _ \ '__|  _  // _ \/ __|
 | |    | | | | | (_| |_| |_| | | | ||  __/ |  | | \ \  __/\__ \
 |_|    |_|_| |_|\__,_|_____|_| |_|\__\___|_|  |_|  \_\___||___/.py
                                                                                                 
     Version: {} | Designed by: {} ({})
 ==============================================================
    """.format(VERSION, AUTHOR, URL))

def get_interacting_residues(pdb_file):
    """
    Core Logic: Analyzes PDB and returns formatted list string or error message.
    """
    try:
        u = mda.Universe(pdb_file)
    except Exception as e:
        return None, f"Error loading PDB: {e}"

    # 1. AUTO-DETECT LIGAND
    forbidden = ["HOH", "SOL", "WAT", "TIP3", "NA", "CL", "MG", "K", "ZN", "CA", "FE", "MN", "PO4", "SO4"]
    selection_str = "not protein and not nucleic"
    for res in forbidden:
        selection_str += f" and not resname {res}"
        
    candidates = u.select_atoms(selection_str)
    
    if len(candidates) == 0:
        return None, "No interacting residues found (No ligand detected)"

    # Heuristic: Ligand is the candidate with the most atoms
    unique_candidates = list(set(candidates.residues))
    best_ligand = None
    max_atoms = 0
    
    for res in unique_candidates:
        n_atoms = len(res.atoms)
        if n_atoms > max_atoms and n_atoms >= 5: 
            max_atoms = n_atoms
            best_ligand = res
            
    if best_ligand is None:
        return None, "No interacting residues found (Only small ions detected)"

    # 2. FIND INTERACTIONS (6.0 Angstroms)
    cutoff = 6.0
    nearby_atoms = u.select_atoms(f"protein and around {cutoff} group lig", lig=best_ligand.atoms)
    unique_residues = list(set(nearby_atoms.residues))
    
    # 3. FORMAT OUTPUT
    unique_residues.sort(key=lambda x: (x.segid, x.resid))
    
    output_list = []
    for res in unique_residues:
        chain_id = str(res.segid) if res.segid else "A"
        res_num = int(res.resid)
        output_list.append((chain_id, res_num))
        
    if not output_list:
        return None, "No interacting residues found"
    
    return str(output_list), None

def save_to_file(content, input_path, output_dir):
    """Helper to save content to <input_name>_InterRes.txt"""
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except Exception as e:
            return False, f"Could not create directory: {e}"

    input_filename = Path(input_path).stem
    output_filename = f"{input_filename}_InterRes.txt"
    full_save_path = os.path.join(output_dir, output_filename)
    
    try:
        with open(full_save_path, "w") as f:
            f.write(content)
        return True, full_save_path
    except Exception as e:
        return False, str(e)

# --- MODES ---

def run_interactive():
    print_header()
    
    # Input
    while True:
        raw_path = input("Please enter the path to your PDB complex file: ").strip()
        pdb_path = raw_path.strip('"').strip("'") # Clean quotes
        if os.path.exists(pdb_path):
            break
        print(f"Error: File not found at: {pdb_path}\n")

    print(f"\nProcessing: {os.path.basename(pdb_path)}...")
    result, error = get_interacting_residues(pdb_path)
    
    if error:
        print(f"\n[!] {error}")
    else:
        print("\n--- Results ---")
        print(result)
        print("---------------\n")
        
        # Save?
        if input("Do you wish to save this output to a text file? (y/n): ").lower().strip() == 'y':
            while True:
                out_dir_raw = input("Enter the directory path to save the file: ").strip()
                out_dir = out_dir_raw.strip('"').strip("'")
                # Attempt save
                success, msg = save_to_file(result, pdb_path, out_dir)
                if success:
                    print(f"\nSuccess! Output saved to:\n{msg}")
                    break
                else:
                    print(f"Error: {msg}. Try again.")
        else:
            print("Output not saved.")
            
    input("\nPress Enter to exit...")

def run_pipeline(args):
    # In pipeline mode, we usually skip the fancy header unless verbose
    if not args.quiet:
        print_header()

    pdb_path = args.input.strip('"').strip("'")
    
    if not os.path.exists(pdb_path):
        print(f"Error: Input file not found: {pdb_path}", file=sys.stderr)
        sys.exit(1)

    result, error = get_interacting_residues(pdb_path)
    
    if error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
        
    # Print result to stdout (for piping)
    print(result)
    
    # Save to file if requested
    if args.output:
        out_dir = args.output.strip('"').strip("'")
        success, msg = save_to_file(result, pdb_path, out_dir)
        if success:
            if not args.quiet: print(f"File saved: {msg}", file=sys.stderr)
        else:
            print(f"Error saving file: {msg}", file=sys.stderr)
            sys.exit(1)

def main():
    # If no arguments are passed, default to Interactive Mode
    if len(sys.argv) == 1:
        run_interactive()
    else:
        # Parse arguments for Pipeline Mode
        parser = argparse.ArgumentParser(description=f"{SCRIPT_NAME} v{VERSION} - Ligand Interaction Finder")
        
        parser.add_argument("-i", "--input", required=True, help="Path to the input PDB complex file")
        parser.add_argument("-o", "--output", help="Directory path to save the output text file (optional)")
        parser.add_argument("-q", "--quiet", action="store_true", help="Suppress header and extra logs (output only results)")
        
        args = parser.parse_args()
        run_pipeline(args)

if __name__ == "__main__":
    main()