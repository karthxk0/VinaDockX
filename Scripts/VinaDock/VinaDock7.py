#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import glob
import time
import datetime
import shutil

# ANSI Color Codes for aesthetic terminal output
class C:
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    MAGENTA = '\033[95m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

def print_banner():
    # Using a raw string (r"") avoids any Invalid Escape Sequence errors
    banner = C.CYAN + r"""
======================================================================                                                       
 __      ___             _____             _    
 \ \    / (_)           |  __ \           | |   
  \ \  / / _ _ __   __ _| |  | | ___   ___| | __
   \ \/ / | | '_ \ / _` | |  | |/ _ \ / __| |/ /
    \  /  | | | | | (_| | |__| | (_) | (__|   < 
     \/   |_|_| |_|\__,_|_____/ \___/ \___|_|\_\.py                                                       
                                                                  
  > Version 2.7  |  > Built by karthxk (https://karthxk0.github.io/)         
  
======================================================================""" + C.RESET
    print(banner)

def check_dependencies():
    required_packages = ["numpy", "scipy", "rdkit", "meeko", "gemmi"]
    missing = []
    for pkg in required_packages:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    return missing

def clean_path(path_str):
    if not path_str:
        return ""
    return path_str.strip().strip('"').strip("'")

def get_vina_executable():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    if sys.platform.startswith('win'):
        vina_path = os.path.join(base_dir, 'vina-win', 'vina_1.2.7_win .exe')
        os_type = "Windows"
    else:
        vina_path = os.path.join(base_dir, 'vina-linux', 'vina_1.2.7_linux_x86_64')
        os_type = "Linux"
        if os.path.exists(vina_path):
            os.chmod(vina_path, 0o755)
            
    return vina_path, os_type

def get_vina_split_executable():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    if sys.platform.startswith('win'):
        search_pattern = os.path.join(base_dir, 'vina_split-win', 'vina_split*.exe')
        matches = glob.glob(search_pattern)
        if matches:
            return matches[0]
        return os.path.join(base_dir, 'vina_split-win', 'vina_split_1.2.7_win.exe')
    else:
        search_pattern = os.path.join(base_dir, 'vina_split-linux', 'vina_1.2.7_linux_x86_64')
        matches = glob.glob(search_pattern)
        split_path = matches[0] if matches else os.path.join(base_dir, 'vina_split-linux', 'vina_split')
        if os.path.exists(split_path):
            os.chmod(split_path, 0o755)
        return split_path

def get_ligands_from_input(input_str):
    ligands = []
    paths = [p.strip() for p in input_str.split(',')]
    for p in paths:
        clean_p = clean_path(p)
        if os.path.isfile(clean_p) and clean_p.lower().endswith('.pdbqt'):
            ligands.append(clean_p)
        elif os.path.isdir(clean_p):
            for root, _, files in os.walk(clean_p):
                for file in files:
                    if file.lower().endswith('.pdbqt'):
                        ligands.append(os.path.join(root, file))
    return list(set(ligands))

def get_search_space(config_path):
    full_space = "Not Specified"
    short_space = ""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            content = f.read().lower()
            if "blind (whole protein)" in content or "blind" in content:
                full_space = "Blind (whole protein)"
                short_space = "Blind"
            elif "targeted (based on residues provided)" in content or "targeted" in content:
                full_space = "Targeted (based on residues provided)"
                short_space = "Targeted"
    except Exception as e:
        pass
    return full_space, short_space

def clean_pdbqt_output(filepath):
    """Removes null bytes (NUL) from the PDBQT output to fix PyMOL multi-model parsing issues."""
    if not os.path.exists(filepath):
        return
    try:
        # Read the file in binary mode to safely capture and remove null bytes
        with open(filepath, 'rb') as f:
            content = f.read()
        
        # Strip out all \x00 (NUL) characters
        cleaned_content = content.replace(b'\x00', b'')
        
        # Write the cleaned binary content back to the file
        with open(filepath, 'wb') as f:
            f.write(cleaned_content)
    except Exception as e:
        print(f"{C.RED}[!] Error scrubbing NUL bytes from PDBQT: {e}{C.RESET}")

def run_docking(vina_path, config_path, receptor_path, flex_path, ligand_path, out_path, log_path, mode_str, search_space_full):
    cmd = [
        vina_path,
        "--config", config_path,
        "--receptor", receptor_path,
        "--out", out_path
    ]
    if flex_path:
        cmd.extend(["--flex", flex_path])
        
    cmd.append("--ligand")
    if isinstance(ligand_path, list):
        cmd.extend(ligand_path)
    else:
        cmd.append(ligand_path)

    header = (
        "# Generated by VinaDock.py v2.8\n"
        "# Built by karthxk (https://karthxk0.github.io/)\n\n"
        f"# Search Space: {search_space_full}\n"
        f"# Mode: {mode_str}\n\n"
        "# Docking Log:\n\n"
    )
    
    success = False
    error_msg = ""
    
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=0
        )
        
        with open(log_path, "wb") as log_file:
            log_file.write(header.encode('utf-8'))
            
            while True:
                char = process.stdout.read(1)
                if not char and process.poll() is not None:
                    break
                if char:
                    sys.stdout.buffer.write(char)
                    sys.stdout.buffer.flush()
                    log_file.write(char)
                    log_file.flush()
            
            if process.returncode == 0:
                success = True
                clean_pdbqt_output(out_path)  # Scrub the NUL bytes automatically
            else:
                error_msg = f"Vina exited with code {process.returncode}"
                
        print()
                
    except Exception as e:
        error_msg = str(e)
        
    return success, error_msg

def generate_sdf_with_meeko(pdbqt_path, sdf_path):
    if sys.platform.startswith('win'):
        meeko_cmd = "mk_export"
    else:
        meeko_cmd = "mk_export.py"
        
    cmd = [meeko_cmd, pdbqt_path, "-s", sdf_path]
    try:
        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if process.returncode == 0:
            return True, ""
        else:
            return False, process.stderr
    except Exception as e:
        return False, str(e)

def split_models(split_exe, pdbqt_file, out_dir_pdbqt, out_dir_sdf, receptor_base, lig_base, space_tag, flex_tag):
    # Pass cwd=out_dir_pdbqt to force vina_split to output in the PDBQT folder instead of the terminal directory
    cmd = [split_exe, "--input", pdbqt_file]
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=out_dir_pdbqt)
    
    base_name = os.path.splitext(os.path.basename(pdbqt_file))[0]
    # Look for the specific _ligand_ outputs generated by Vina Split
    pattern = os.path.join(out_dir_pdbqt, f"{base_name}_ligand_*.pdbqt")
    split_files = glob.glob(pattern)
    
    split_pdbqt_dir = os.path.join(out_dir_pdbqt, "Split Models")
    split_sdf_dir = os.path.join(out_dir_sdf, "Split Models")
    
    sp_pdbqt_count = 0
    sp_sdf_count = 0
    
    for sf in split_files:
        # Extract numbering (e.g. "01" from output_ligand_01.pdbqt)
        part = sf.split('_ligand_')[-1].replace('.pdbqt', '')
        try:
            model_num = int(part) # Strips leading zeroes
        except ValueError:
            model_num = part
            
        new_name_base = f"{receptor_base}-{lig_base}-Model{model_num}{space_tag}{flex_tag}_VinaDock_Output"
        new_pdbqt = os.path.join(split_pdbqt_dir, f"{new_name_base}.pdbqt")
        new_sdf = os.path.join(split_sdf_dir, f"{new_name_base}.sdf")
        
        shutil.move(sf, new_pdbqt)
        sp_pdbqt_count += 1
        
        # Move flexible residues if they exist
        flex_file = os.path.join(out_dir_pdbqt, f"{base_name}_flex_{part}.pdbqt")
        if os.path.exists(flex_file):
            new_flex_name = f"{receptor_base}-{lig_base}-Model{model_num}{space_tag}-FlexRes_VinaDock_Output.pdbqt"
            shutil.move(flex_file, os.path.join(split_pdbqt_dir, new_flex_name))
        
        succ, err = generate_sdf_with_meeko(new_pdbqt, new_sdf)
        if succ:
            sp_sdf_count += 1
            
    return sp_pdbqt_count, sp_sdf_count

def generate_merged_log(log_files, output_dir, receptor_name, merge_filename, mode_str, search_space_full):
    if not log_files:
        return
    
    merged_path = os.path.join(output_dir, merge_filename)
    
    with open(merged_path, 'w') as out_f:
        out_f.write("# Generated by VinaDock.py v2.8\n")
        out_f.write("# Built by karthxk (https://karthxk0.github.io/)\n\n")
        out_f.write(f"# Search Space: {search_space_full}\n")
        out_f.write(f"# Mode: {mode_str}\n\n")
        
        first_log = log_files[0]
        try:
            with open(first_log, 'r') as f:
                first_lines = f.readlines()
                
            vina_start = 0
            verbosity_end = 0
            
            for i, line in enumerate(first_lines):
                if line.startswith("AutoDock Vina"):
                    vina_start = i
                if line.startswith("Verbosity:"):
                    verbosity_end = i
                    break
            
            if vina_start > 0 and verbosity_end > 0:
                for i in range(vina_start, verbosity_end + 1):
                    if not first_lines[i].startswith("Ligand:"):
                        out_f.write(first_lines[i])
                out_f.write("\n=================================================================\n\n")
        except Exception as e:
            print(f"{C.RED}[!] Error reading baseline log for merging: {e}{C.RESET}")
            return

        for log_path in log_files:
            try:
                with open(log_path, 'r') as f:
                    lines = f.readlines()
                    
                ligand_line = ""
                seed_line = ""
                table_start = -1
                
                for i, line in enumerate(lines):
                    if line.startswith("Ligand:"):
                        ligand_line = line
                    elif "Performing docking (random seed:" in line:
                        seed = line.split("random seed:")[1].split(")")[0].strip()
                        seed_line = f"random seed: {seed}\n\n"
                    elif line.startswith("mode |"):
                        table_start = i
                        
                if table_start != -1:
                    if ligand_line: out_f.write(ligand_line + "\n")
                    if seed_line: out_f.write(seed_line)
                    
                    for i in range(table_start, len(lines)):
                        out_f.write(lines[i])
                        
                    out_f.write("\n=================================================================\n\n")
            except Exception as e:
                print(f"{C.RED}[!] Error extracting data from {log_path}: {e}{C.RESET}")

def interactive_mode():
    print(f"{C.GREEN}[?] Select Docking Mode:{C.RESET}")
    print("1. Standard Docking")
    print("2. Flexible Docking")
    mode_choice = input(f"{C.CYAN}Enter 1 or 2: {C.RESET}").strip()
    
    is_flex = (mode_choice == '2')
    
    if is_flex:
        receptor = clean_path(input(f"{C.CYAN}Enter path to Prepared Receptor (Rigid PDBQT): {C.RESET}"))
        flex_res = clean_path(input(f"{C.CYAN}Enter path to Flexible Residues (Flex PDBQT): {C.RESET}"))
    else:
        receptor = clean_path(input(f"{C.CYAN}Enter path to Prepared Receptor (PDBQT): {C.RESET}"))
        flex_res = None
        
    ligands_in = input(f"{C.CYAN}Enter path to Prepared Ligands (Files/Folders, comma separated): {C.RESET}")
    ligands = get_ligands_from_input(ligands_in)
    
    simultaneous_docking = False
    if len(ligands) > 1:
        print(f"\n{C.CYAN}Do you wish to perform Simultaneous multiple ligand docking?{C.RESET}")
        print("1. Yes")
        print("2. No")
        sim_choice = input(f"{C.CYAN}Enter 1 or 2: {C.RESET}").strip()
        simultaneous_docking = (sim_choice == '1')
    
    config = clean_path(input(f"{C.CYAN}Enter path to Config (TXT) file: {C.RESET}"))
    out_dir = clean_path(input(f"{C.CYAN}Enter directory to save outputs: {C.RESET}"))
    
    merge_logs = False
    if len(ligands) > 1 and not simultaneous_docking:
        merge_in = input(f"{C.CYAN}Merge all logs into a single file? (1 for Yes, 2 for No): {C.RESET}").strip()
        merge_logs = (merge_in == '1' or merge_in.lower() in ['y', 'yes'])
    
    return receptor, flex_res, ligands, config, out_dir, merge_logs, is_flex, simultaneous_docking

def main():
    print_banner()
    
    print(f"{C.MAGENTA}[*] Checking system dependencies...{C.RESET}")
    missing_pkgs = check_dependencies()
    if missing_pkgs:
        print(f"{C.RED}[!] Critical Error: Missing required Python packages: {', '.join(missing_pkgs)}{C.RESET}")
        print(f"{C.YELLOW}[i] Please install them using: pip install {(' '.join(missing_pkgs))}{C.RESET}")
        sys.exit(1)
    print(f"{C.GREEN}[+] Python dependencies verified.{C.RESET}")
    
    vina_path, os_type = get_vina_executable()
    print(f"{C.GREEN}[+] Detected OS: {os_type}{C.RESET}")
    
    if not os.path.exists(vina_path):
        print(f"{C.RED}[!] Critical Error: Vina executable not found at expected path: {vina_path}{C.RESET}")
        sys.exit(1)
    print(f"{C.GREEN}[+] Vina executable verified.{C.RESET}\n")

    parser = argparse.ArgumentParser(description="VinaDock Automation Script")
    parser.add_argument('--receptor', help='Path to rigid receptor PDBQT')
    parser.add_argument('--flex', help='Path to flexible residues PDBQT')
    parser.add_argument('--ligands', help='Comma-separated paths to ligand files or folders')
    parser.add_argument('--config', help='Path to Config TXT file')
    parser.add_argument('--outdir', help='Output directory')
    parser.add_argument('--merge', action='store_true', help='Merge log files')
    parser.add_argument('--simultaneous', action='store_true', help='Run simultaneous multiple ligand docking')
    parser.add_argument('--cli', action='store_true', help='Run in CLI mode without prompts')
    
    args = parser.parse_args()

    if args.cli or (args.receptor and args.ligands and args.config and args.outdir):
        print(f"{C.CYAN}[*] Running in CLI Arguments Mode...{C.RESET}\n")
        receptor = clean_path(args.receptor)
        flex_res = clean_path(args.flex) if args.flex else None
        ligands = get_ligands_from_input(args.ligands)
        config = clean_path(args.config)
        out_dir = clean_path(args.outdir)
        merge_logs = args.merge
        simultaneous_docking = args.simultaneous
        is_flex = bool(flex_res)
    else:
        receptor, flex_res, ligands, config, out_dir, merge_logs, is_flex, simultaneous_docking = interactive_mode()

    if not os.path.exists(receptor):
        print(f"{C.RED}[!] Receptor file not found: {receptor}{C.RESET}")
        sys.exit(1)
    if not os.path.exists(config):
        print(f"{C.RED}[!] Config file not found: {config}{C.RESET}")
        sys.exit(1)
        
    if os.path.isfile(out_dir):
        print(f"{C.RED}[!] Error: You provided a file path for the output directory ({os.path.basename(out_dir)}).{C.RESET}")
        print(f"{C.YELLOW}[i] Please run the script again and provide a folder path instead.{C.RESET}")
        sys.exit(1)
        
    os.makedirs(out_dir, exist_ok=True)
    
    # Check for Vina Split executable
    split_exe = get_vina_split_executable()
    can_split = os.path.exists(split_exe)
    if not can_split:
        print(f"{C.YELLOW}[i] vina_split executable not found. Model splitting will be skipped.{C.RESET}\n")
    
    # Create required subdirectories
    pdbqt_dir = os.path.join(out_dir, "PDBQT Output")
    sdf_dir = os.path.join(out_dir, "SDF Output")
    log_dir = os.path.join(out_dir, "Log")
    
    os.makedirs(pdbqt_dir, exist_ok=True)
    os.makedirs(sdf_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    if can_split:
        os.makedirs(os.path.join(pdbqt_dir, "Split Models"), exist_ok=True)
        os.makedirs(os.path.join(sdf_dir, "Split Models"), exist_ok=True)
    
    if not ligands:
        print(f"{C.RED}[!] No valid ligand PDBQT files found from input.{C.RESET}")
        sys.exit(1)

    mode_str = "Flexible Docking" if is_flex else "Standard Docking"
    receptor_base = os.path.splitext(os.path.basename(receptor))[0]
    
    search_space_full, search_space_short = get_search_space(config)
    space_tag = f"-{search_space_short}" if search_space_short else ""
    
    print(f"\n{C.BOLD}{C.MAGENTA}================= DOCKING INITIALIZED ================={C.RESET}")
    print(f"  {C.CYAN}Mode:{C.RESET} {mode_str}")
    print(f"  {C.CYAN}Search Space:{C.RESET} {search_space_full}")
    print(f"  {C.CYAN}Receptor:{C.RESET} {receptor_base}")
    if simultaneous_docking and len(ligands) > 1:
        print(f"  {C.CYAN}Ligands to process:{C.RESET} {len(ligands)} (Simultaneous Docking)")
    else:
        print(f"  {C.CYAN}Ligands to process:{C.RESET} {len(ligands)}")
    print(f"  {C.CYAN}Output Directory:{C.RESET} {out_dir}")
    print(f"{C.BOLD}{C.MAGENTA}======================================================={C.RESET}\n")

    success_count = 0
    sdf_count = 0
    total_split_pdbqt = 0
    total_split_sdf = 0
    fail_count = 0
    failed_details = []
    generated_logs = []
    
    # Determine processing batches based on simultaneous parameter
    if simultaneous_docking and len(ligands) > 1:
        batches = [(ligands, "MultiLigand")]
    else:
        batches = [([lig], os.path.splitext(os.path.basename(lig))[0]) for lig in ligands]
    
    for idx, (batch_ligands, lig_base) in enumerate(batches, 1):
        flex_tag = "-Flex" if is_flex else ""
        
        out_filename = f"{receptor_base}-{lig_base}{space_tag}{flex_tag}_VinaDock_Output.pdbqt"
        log_filename = f"{receptor_base}-{lig_base}{space_tag}{flex_tag}_VinaDock_Log.txt"
        
        # Route files to their respective subdirectories
        out_path = os.path.join(pdbqt_dir, out_filename)
        log_path = os.path.join(log_dir, log_filename)
        
        if simultaneous_docking and len(ligands) > 1:
            print(f"{C.GREEN}[> Processing Simultaneous Batch] {C.CYAN}Total Ligands: {len(batch_ligands)}{C.RESET}")
        else:
            print(f"{C.GREEN}[> Processing {idx}/{len(batches)}] {C.CYAN}Ligand: {lig_base}{C.RESET}")
        
        success, error = run_docking(vina_path, config, receptor, flex_res, batch_ligands, out_path, log_path, mode_str, search_space_full)
        
        if success:
            print(f"{C.GREEN}[+] Successfully docked: {lig_base}{C.RESET}")
            success_count += 1
            generated_logs.append(log_path)
            
            # --- Meeko SDF Generation Step (Combined File) ---
            print(f"{C.YELLOW}[*] Converting combined output to SDF using Meeko...{C.RESET}")
            sdf_filename = out_filename.replace(".pdbqt", ".sdf")
            sdf_path = os.path.join(sdf_dir, sdf_filename)
            sdf_success, sdf_error = generate_sdf_with_meeko(out_path, sdf_path)
            
            if sdf_success:
                print(f"{C.GREEN}[+] Successfully generated SDF: {sdf_filename}{C.RESET}")
                sdf_count += 1
            else:
                print(f"{C.RED}[!] Failed to generate SDF for {lig_base}: {sdf_error}{C.RESET}")
                
            # --- Vina Split Step ---
            if can_split:
                print(f"{C.YELLOW}[*] Splitting PDBQT output into multiple individual models...{C.RESET}")
                sp_pdbqt, sp_sdf = split_models(split_exe, out_path, pdbqt_dir, sdf_dir, receptor_base, lig_base, space_tag, flex_tag)
                total_split_pdbqt += sp_pdbqt
                total_split_sdf += sp_sdf
                if sp_pdbqt > 0:
                    print(f"{C.GREEN}[+] Successfully generated {sp_pdbqt} split PDBQTs and {sp_sdf} split SDFs.{C.RESET}\n")
                else:
                    print(f"{C.RED}[!] Splitting failed or produced no output models.{C.RESET}\n")
            else:
                print("") # Spacing
                
        else:
            print(f"{C.RED}[!] Failed to dock: {receptor_base} + {lig_base}{C.RESET}")
            print(f"{C.RED}    Error: {error}{C.RESET}\n")
            fail_count += 1
            failed_details.append(f"{receptor_base} + {lig_base} ({error})")

    merged_log_name = ""
    if merge_logs and len(generated_logs) > 1 and not simultaneous_docking:
        timestamp = datetime.datetime.now().strftime("%d%m%Y-%H%M")
        merged_log_name = f"MergedLogs-{receptor_base}{space_tag}{flex_tag}_VinaDock_{timestamp}.txt"
        print(f"{C.CYAN}[*] Merging logs into: {merged_log_name}...{C.RESET}")
        generate_merged_log(generated_logs, log_dir, receptor_base, merged_log_name, mode_str, search_space_full)

    print(f"\n{C.BOLD}{C.MAGENTA}==================== FINAL SUMMARY ===================={C.RESET}")
    print(f"{C.CYAN}Mode:{C.RESET} {mode_str}")
    print(f"{C.CYAN}Search Space:{C.RESET} {search_space_full}")
    print(f"{C.CYAN}Total Inputs:{C.RESET} {len(ligands)}")
    print(f"{C.GREEN}Successful Dockings:{C.RESET} {success_count}")
    print(f"{C.RED}Failed Dockings:{C.RESET} {fail_count}")
    
    if fail_count > 0:
        print(f"\n{C.RED}Failures:{C.RESET}")
        for f in failed_details:
            print(f"  - {f}")
            
    print(f"\n{C.CYAN}Files Generated:{C.RESET}")
    print(f"  - Output PDBQTs (Combined): {success_count}")
    print(f"  - Output SDFs (Combined): {sdf_count}")
    if can_split:
        print(f"  - Output PDBQTs (Split): {total_split_pdbqt}")
        print(f"  - Output SDFs (Split): {total_split_sdf}")
    print(f"  - Individual Logs: {success_count}")
    if merged_log_name:
        print(f"  - Merged Log: 1 ({merged_log_name})")
        
    print(f"\n{C.GREEN}All files have been saved to subdirectories within:{C.RESET} {os.path.abspath(out_dir)}")
    print(f"{C.BOLD}{C.MAGENTA}======================================================={C.RESET}\n")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n{C.RED}[!] Process interrupted by user. Exiting...{C.RESET}")
        sys.exit(1)