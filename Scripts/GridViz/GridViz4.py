import sys
import os
import argparse

# --- HELPER FUNCTIONS ---

def print_banner():
    """Prints the aesthetic header."""
    print(r"""
   
=========================================================================                                                               
   _____      _     ___      ___     
  / ____|    (_)   | \ \    / (_)    
 | |  __ _ __ _  __| |\ \  / / _ ____
 | | |_ | '__| |/ _` | \ \/ / | |_  /
 | |__| | |  | | (_| |  \  /  | |/ / 
  \_____|_|  |_|\__,_|   \/   |_/___|.py  
                                                                  
  > Version 2.4  |  > Designed by karthxk (https://karthxk0.github.io/)         
  
=========================================================================  
    """)
    print("\n")

def get_clean_path(prompt_text):
    """Asks for a file path, strips quotes, and normalizes slashes."""
    raw_input = input(prompt_text).strip()
    if (raw_input.startswith('"') and raw_input.endswith('"')) or \
       (raw_input.startswith("'") and raw_input.endswith("'")):
        clean_input = raw_input[1:-1]
    else:
        clean_input = raw_input
    return os.path.normpath(clean_input)

def parse_list_string(input_str):
    """Parses '[10.5, 20.0, 30.0]' into a list of floats."""
    clean_str = input_str.replace('[', '').replace(']', '').strip()
    try:
        parts = [float(x.strip()) for x in clean_str.split(',')]
        if len(parts) != 3: raise ValueError
        return parts
    except ValueError:
        return None

def parse_autodock_config(config_path):
    """Reads an AutoDock config file and extracts center and size."""
    data = {}
    try:
        with open(config_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or "=" not in line: continue
                parts = line.split('=')
                key = parts[0].strip()
                try:
                    value = float(parts[1].strip())
                    data[key] = value
                except ValueError: continue

        required = ['center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z']
        if all(k in data for k in required):
            center = [data['center_x'], data['center_y'], data['center_z']]
            size = [data['size_x'], data['size_y'], data['size_z']]
            return center, size
        else:
            return None, None
    except FileNotFoundError:
        return None, None

# --- CORE LOGIC ---

def process_pdb(input_pdb, output_dir, center, size, silent=False):
    """Generates the Grid PDB and the PML launcher."""
    
    if not os.path.exists(input_pdb):
        if not silent: print(f"  ! Error: Input PDB file not found: {input_pdb}")
        sys.exit(1)

    # Prepare Filenames
    base_name = os.path.basename(input_pdb)
    root_name = os.path.splitext(base_name)[0]
    out_pdb_name = f"{root_name}_grid.pdb"
    out_pml_name = f"{root_name}_view.pml"
    final_pdb_path = os.path.join(output_dir, out_pdb_name)
    final_pml_path = os.path.join(output_dir, out_pml_name)

    # Create Box Geometry
    cx, cy, cz = center
    sx, sy, sz = size
    hx, hy, hz = sx/2.0, sy/2.0, sz/2.0
    
    vertices = [
        [cx-hx, cy-hy, cz-hz], [cx+hx, cy-hy, cz-hz],
        [cx+hx, cy+hy, cz-hz], [cx-hx, cy+hy, cz-hz],
        [cx-hx, cy-hy, cz+hz], [cx+hx, cy-hy, cz+hz],
        [cx+hx, cy+hy, cz+hz], [cx-hx, cy+hy, cz+hz],
    ]
    connections = [
        (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7)
    ]

    # Read PDB
    lines = []
    last_serial = 0
    with open(input_pdb, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    s = int(line[6:11])
                    if s > last_serial: last_serial = s
                except ValueError: pass
            if not line.startswith("END"):
                lines.append(line)

    # Create Box Atoms/Bonds
    box_lines = []
    start_serial = last_serial + 1
    atom_serials = []
    for i, (vx, vy, vz) in enumerate(vertices):
        serial = start_serial + i
        atom_serials.append(serial)
        box_lines.append(f"HETATM{serial:5d}  C   BOX Z   1    {vx:8.3f}{vy:8.3f}{vz:8.3f}  1.00  0.00           C  \n")

    conect_lines = []
    for s_idx, e_idx in connections:
        conect_lines.append(f"CONECT{atom_serials[s_idx]:5d}{atom_serials[e_idx]:5d}\n")

    # Write PDB
    try:
        with open(final_pdb_path, 'w') as f:
            f.writelines(lines)
            f.write(f"\nREMARK 999 GRID BOX: Center={center} Size={size}\n")
            f.writelines(box_lines)
            f.writelines(conect_lines)
            f.write("END\n")
    except Exception as e:
        if not silent: print(f"  ! Error writing PDB: {e}")
        sys.exit(1)

    # Write PML
    pml_content = f"""load {out_pdb_name}
hide everything
show cartoon, polymer
show sticks, hetatm
color blue, polymer
color red, hetatm
color white, resn BOX
set stick_radius, 0.15, resn BOX
zoom
"""
    try:
        with open(final_pml_path, 'w') as f:
            f.write(pml_content)
    except Exception as e:
        if not silent: print(f"  ! Error writing PML: {e}")
        sys.exit(1)

    if not silent:
        print("\n" + "="*60)
        print(f"  SUCCESS!")
        print(f"  Files saved to: {output_dir}")
        print(f"  1. PDB File: {out_pdb_name}")
        print(f"  2. Launcher: {out_pml_name} (Open this file in PyMOL)")
        print("="*60 + "\n")

# --- MODES ---

def run_interactive():
    print_banner()
    
    # 1. PDB Input
    pdb_path = get_clean_path("  [?] Enter path to input PDB file: ")
    while not os.path.exists(pdb_path):
        print("      ! File not found.")
        pdb_path = get_clean_path("  [?] Enter path to input PDB file: ")

    # 2. Grid Input
    print("\n  [?] How would you like to input grid dimensions?")
    print("      1. Direct Input (Enter center and size manually)")
    print("      2. AutoDock Config File (Read from .txt/.gpf)")
    choice = input("      > Select (1/2): ").strip()
    
    center, size = None, None
    if choice == '1':
        while center is None:
            c_str = input("  [?] Enter Grid CENTER [x, y, z]: ")
            center = parse_list_string(c_str)
            if center is None: print("      ! Invalid format. Use [x, y, z]")
        while size is None:
            s_str = input("  [?] Enter Grid SIZE   [x, y, z]: ")
            size = parse_list_string(s_str)
            if size is None: print("      ! Invalid format. Use [x, y, z]")
    elif choice == '2':
        while center is None:
            cfg_path = get_clean_path("  [?] Enter path to Config file: ")
            center, size = parse_autodock_config(cfg_path)
            if center is None: print("      Let's try again.")
    else:
        print("  ! Invalid choice.")
        return

    # 3. Output Dir
    out_dir = get_clean_path("\n  [?] Enter directory to save output files: ")
    if not os.path.exists(out_dir):
        create = input(f"      Directory '{out_dir}' does not exist. Create it? (y/n): ")
        if create.lower().startswith('y'):
            os.makedirs(out_dir)
        else:
            return

    print(f"\n  ... Processing ...")
    process_pdb(pdb_path, out_dir, center, size)

def run_cli():
    parser = argparse.ArgumentParser(description="Gridmaker: Embed AutoDock grids into PDB files.")
    
    # Required input
    parser.add_argument("-p", "--pdb", required=True, help="Input PDB file path")
    parser.add_argument("-o", "--out", required=True, help="Output directory path")
    
    # Grid definition (Mutually exclusive group)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--config", help="AutoDock config file path")
    group.add_argument("-m", "--manual", action="store_true", help="Flag to indicate manual center/size input")

    # Manual arguments (only needed if -m is used)
    # We accept strings now, e.g. "[-10, 20, 30]"
    parser.add_argument("--center", type=str, help="Grid Center e.g. '[-55.9, 5.4, 72.2]'")
    parser.add_argument("--size", type=str, help="Grid Size e.g. '[18.3, 21.4, 20.9]'")

    args = parser.parse_args()

    # Logic
    center, size = None, None

    if args.config:
        center, size = parse_autodock_config(args.config)
        if center is None:
            print("Error: Could not parse config file.")
            sys.exit(1)
    
    elif args.manual:
        if not args.center or not args.size:
            print("Error: In manual mode (-m), you must provide --center and --size")
            sys.exit(1)
            
        center = parse_list_string(args.center)
        size = parse_list_string(args.size)
        
        if center is None or size is None:
            print("Error: Invalid format for center or size. Use '[x, y, z]' including brackets.")
            sys.exit(1)

    # Output Dir handling
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    process_pdb(args.pdb, args.out, center, size, silent=True)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_cli()
    else:
        run_interactive()