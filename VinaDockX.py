import os
import sys
import argparse
import subprocess
import datetime
import shutil
import glob
import re
import time
import importlib.util

# ==========================================
# Dependency Management & Aesthetics
# ==========================================
def install_and_import_aesthetics():
    """Ensures colorama and tqdm are available for the master script."""
    missing = []
    for pkg in ['colorama', 'tqdm']:
        if importlib.util.find_spec(pkg) is None:
            missing.append(pkg)
    if missing:
        print(f"\n[*] Installing required master dependencies: {', '.join(missing)}...")
        subprocess.run([sys.executable, "-m", "pip", "install", *missing], stdout=subprocess.DEVNULL)
        
install_and_import_aesthetics()
from colorama import init, Fore, Style
from tqdm import tqdm
init(autoreset=True)

class C:
    CYAN = Fore.LIGHTCYAN_EX
    GREEN = Fore.LIGHTGREEN_EX
    YELLOW = Fore.LIGHTYELLOW_EX
    RED = Fore.LIGHTRED_EX
    MAGENTA = Fore.LIGHTMAGENTA_EX
    WHITE = Fore.LIGHTWHITE_EX
    RESET = Style.RESET_ALL
    BOLD = Style.BRIGHT

def print_banner(log_file=None, to_terminal=True):
    banner = f"""{C.CYAN}
=============================================================================================  

░██    ░██ ░██                      ░███████                         ░██       ░██    ░██ 
░██    ░██                          ░██   ░██                        ░██        ░██  ░██  
░██    ░██ ░██░████████   ░██████   ░██    ░██  ░███████   ░███████  ░██    ░██  ░██░██   
░██    ░██ ░██░██    ░██       ░██  ░██    ░██ ░██    ░██ ░██    ░██ ░██   ░██    ░███    
 ░██  ░██  ░██░██    ░██  ░███████  ░██    ░██ ░██    ░██ ░██        ░███████    ░██░██   
  ░██░██   ░██░██    ░██ ░██   ░██  ░██   ░██  ░██    ░██ ░██    ░██ ░██   ░██  ░██  ░██  
   ░███    ░██░██    ░██  ░█████░██ ░███████    ░███████   ░███████  ░██    ░██░██    ░██ 
                                                                                                                                                                         
  > Version 2.9  |  > Designed by karthxk (https://karthxk0.github.io/)          
  
============================================================================================={C.RESET}
"""
    if to_terminal:
        print(banner)
    if log_file:
        clean_banner = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])').sub('', banner)
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(clean_banner + "\n")

# ==========================================
# Helper Functions
# ==========================================
def clean_input(prompt_text, sample=""):
    print(f"\n{C.YELLOW}{prompt_text}{C.RESET}")
    if sample:
        print(f"{C.CYAN}   Example: {sample}{C.RESET}")
    return input(f"{C.GREEN}> {C.RESET}").strip().strip("\"'")

def log_msg(msg, log_file=None, use_tqdm=False):
    if use_tqdm:
        tqdm.write(msg)
    else:
        print(msg)
        
    if log_file:
        clean_msg = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])').sub('', msg)
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(clean_msg + "\n")

def check_pipeline_dependencies(steps):
    print(f"{C.MAGENTA}[*] Verifying Sub-Script Dependencies...{C.RESET}")
    reqs = set()
    if 1 in steps: reqs.add('MDAnalysis')
    if 2 in steps: reqs.update(['numpy', 'Bio']) 
    if 4 in steps: reqs.update(['meeko', 'tqdm', 'Bio', 'numpy', 'scipy', 'rdkit', 'gemmi'])
    if 5 in steps: reqs.update(['rdkit', 'tqdm', 'numpy', 'scipy', 'rdkit', 'meeko', 'gemmi'])
    # Updated dependencies based on VinaDock v2.6 additions
    if 6 in steps: reqs.update(['numpy', 'scipy', 'rdkit', 'meeko', 'gemmi'])
    if 7 in steps: reqs.update(['pandas', 'openpyxl', 'colorama', 'tqdm'])
    
    missing = []
    for pkg in reqs:
        if importlib.util.find_spec(pkg) is None:
            missing.append('biopython' if pkg == 'Bio' else pkg)
            
    if missing:
        print(f"{C.RED}[!] Missing required packages for selected pipeline steps: {', '.join(missing)}{C.RESET}")
        print(f"{C.RED}    Please run: pip install {' '.join(missing)}{C.RESET}")
        sys.exit(1)
    print(f"{C.GREEN}[+] All dependencies satisfied.{C.RESET}")

def extract_residues(input_val):
    if not input_val: return None
    content = input_val
    if os.path.isfile(input_val):
        try:
            with open(input_val, 'r') as f: content = f.read()
        except Exception: pass

    residues = []
    for chain, resnum in re.findall(r'([A-Za-z])\s*:\s*(\d+)', content):
        residues.append(f"{chain.upper()}:{resnum}")
    for chain, resnum in re.findall(r"\(\s*['\"]([A-Za-z])['\"]\s*,\s*(\d+)\s*\)", content):
        residues.append(f"{chain.upper()}:{resnum}")
    
    unique_res = list(dict.fromkeys(residues))
    return ",".join(unique_res)

def find_sub_scripts():
    scripts_base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Scripts")
    script_names = ['InterResFi', 'GridConfigGen', 'GridViz', 'PrepProt', 'PrepLig', 'VinaDock', 'BindReSort']
    paths = {}
    for s in script_names:
        py_files = glob.glob(os.path.join(scripts_base, s, "*.py"))
        paths[s] = py_files[0] if py_files else None
    
    missing = [k for k, v in paths.items() if v is None]
    if missing:
        print(f"{C.RED}[!] Missing scripts: {', '.join(missing)}{C.RESET}")
        print(f"{C.RED}    Ensure folder structure is exactly ./Scripts/<Name>/<script>.py{C.RESET}")
        sys.exit(1)
    return paths

def parse_run_sequence(seq_str):
    steps = set()
    parts = seq_str.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            try:
                start, end = map(int, part.split('-'))
                steps.update(range(start, end + 1))
            except ValueError: pass
        elif part.isdigit():
            steps.add(int(part))
    return sorted(list(steps))

# ==========================================
# Input Logic & Sequencing
# ==========================================
def gather_inputs(steps):
    cfg = {}
    print(f"\n{C.BOLD}{C.MAGENTA}--- Configuring Pipeline Sequence ---{C.RESET}")
    
    cfg['base_out'] = clean_input("Enter the path to the directory where you want to save the Output:")
    
    if 2 in steps:
        print(f"\n{C.BOLD}{C.WHITE}[Search Space Selection]{C.RESET}")
        print("  1. Blind (whole protein)")
        print("  2. Targeted (based on specific residues)")
        cfg['search_space'] = 'blind' if clean_input("Select search space (1 or 2):") == '1' else 'targeted'
        
        if cfg['search_space'] == 'blind':
            cfg['receptor_pdb'] = clean_input("Enter Receptor (PDB) file path (Accept .pdb files only):")
            
        elif cfg['search_space'] == 'targeted':
            print(f"\n{C.BOLD}{C.WHITE}[Targeted Residue Input]{C.RESET}")
            print("  1. Autodetect interacting residues (6Å) from co-crystallized structure")
            print("  2. Input residues manually")
            t_choice = clean_input("Select targeted input method (1 or 2):")
            
            if t_choice == '1':
                cfg['targeted_auto'] = True
                cfg['receptor_pdb'] = clean_input("Enter the co-crystallized structure (PDB) file path:")
            else:
                cfg['targeted_auto'] = False
                cfg['receptor_pdb'] = clean_input("Enter Receptor (PDB) file path:")
                cfg['manual_res_raw'] = clean_input("Enter residues list OR path to .txt file:", "[('A', 12), ('B', 15)] OR A:12, B:15")
                cfg['manual_res_clean'] = extract_residues(cfg['manual_res_raw'])

    if 'receptor_pdb' not in cfg and any(s in steps for s in [1, 4]):
        cfg['receptor_pdb'] = clean_input("Enter the main Receptor (PDB) file:")

    if 3 in steps and 2 not in steps:
        cfg['config_txt'] = clean_input("Step 3 selected but Step 2 skipped. Enter path to Config (.txt) file:")

    if 4 in steps:
        print(f"\n{C.BOLD}{C.WHITE}[Protein Preparation]{C.RESET}")
        print("  1. Standard Docking (Rigid)")
        print("  2. Flexible Docking")
        cfg['docking_type'] = 'standard' if clean_input("Select docking type (1 or 2):") == '1' else 'flexible'
        
        if cfg['docking_type'] == 'flexible':
            print(f"\n{C.BOLD}{C.WHITE}[Flexible Residues Input]{C.RESET}")
            print("  1. Input flexible residues manually")
            print("  2. Set interacting residues (6Å) from co-crystallized structure as flexible")
            flex_choice = clean_input("Select flexible input method (1 or 2):")
            
            if flex_choice == '1':
                cfg['flex_res_raw'] = clean_input("Enter residues list OR path to .txt file:", "A:12, B:15")
                cfg['flex_res_clean'] = extract_residues(cfg['flex_res_raw'])
            else:
                cfg['flex_auto'] = True
                if cfg.get('targeted_auto') and (1 in steps):
                    print(f"{C.GREEN}   [*] Will auto-use interacting residues calculated in earlier steps.{C.RESET}")
                else:
                    print(f"\n{C.BOLD}{C.WHITE}[Co-crystallized Structure Source]{C.RESET}")
                    print("  1. From input receptor file (Use if the initial receptor was a complex)")
                    print("  2. Input new co-crystallized structure")
                    co_choice = clean_input("Select source (1 or 2):")
                    if co_choice == '2':
                        cfg['prep_cocrystal'] = clean_input("Enter the new receptor-ligand cocrystallized PDB file:")
                    else:
                        cfg['prep_cocrystal'] = cfg['receptor_pdb']

    if 5 in steps:
        cfg['ligands_input'] = clean_input("Enter ligand file(s) (SDF/MOL2) or folder(s) separated by comma:")
        
    if 6 in steps:
        if 4 not in steps:
            cfg['vina_receptor'] = clean_input("Step 6 selected but Step 4 skipped. Enter Receptor PDBQT file:")
            if cfg.get('docking_type', 'standard') == 'flexible':
                cfg['vina_flex'] = clean_input("Enter Flexible Residues PDBQT file:")
        if 5 not in steps:
            cfg['vina_ligands'] = clean_input("Step 6 selected but Step 5 skipped. Enter Ligands folder/files (.pdbqt):")
        if 2 not in steps and 3 not in steps and 'config_txt' not in cfg:
            cfg['config_txt'] = clean_input("Step 6 selected but Step 2 skipped. Enter Config (.txt) file:")

    if 7 in steps and 6 not in steps:
        cfg['bindresort_input'] = clean_input("Step 7 selected but Step 6 skipped. Enter Vina log (.txt) files/folders:")
        
    if 3 in steps and 'receptor_pdb' not in cfg and 'vina_receptor' not in cfg:
        pdb_opt = clean_input("Step 3 selected. Enter Receptor file (PDB/PDBQT) for Grid Visualization (Press Enter to skip):")
        if pdb_opt:
            cfg['grid_receptor'] = pdb_opt

    # --- INTELLIGENT PIPELINE PRUNING ---
    if 1 in steps:
        needs_interres = False
        if cfg.get('search_space') == 'targeted' and cfg.get('targeted_auto'): needs_interres = True
        if cfg.get('docking_type') == 'flexible' and cfg.get('flex_auto'): needs_interres = True
        
        if not needs_interres and len(steps) > 1:
            steps.remove(1)
            print(f"\n{C.YELLOW}[*] Auto-Optimization: Step 1 (InterResFi) will be skipped as it's not required for Blind/Standard parameters.{C.RESET}")

    target_file = cfg.get('receptor_pdb') or cfg.get('vina_receptor') or cfg.get('grid_receptor') or cfg.get('bindresort_input') or cfg.get('ligands_input') or "Target"
    target_file = target_file.split(',')[0].strip()
    base_name = os.path.splitext(os.path.basename(target_file))[0]
    base_name = base_name.replace('_PrepProt', '').replace('_rigid', '').replace('_flex', '')
    
    ts = datetime.datetime.now().strftime("%d%m%Y-%H%M")
    cfg['main_out_dir'] = os.path.join(cfg['base_out'], f"{base_name}_VinaDockX_{ts}")
    cfg['log_file'] = os.path.join(cfg['main_out_dir'], f"{base_name}_VinaDockX-Log_{ts}.txt")
    cfg['run_id'] = ts

    return cfg

def confirm_inputs(cfg, steps):
    summary = f"\n{C.BOLD}{C.MAGENTA}================================================================={C.RESET}\n"
    summary += f"{C.BOLD}{C.WHITE}                     PIPELINE CONFIGURATION CONFIRMATION{C.RESET}\n"
    summary += f"{C.BOLD}{C.MAGENTA}================================================================={C.RESET}\n"
    
    summary += f"{C.CYAN}Scripts Selected:{C.RESET} {steps}\n"
    summary += f"{C.CYAN}Base Output Dir: {C.RESET} {cfg.get('base_out')}\n"
    summary += f"{C.CYAN}Main Session Dir:{C.RESET} {cfg.get('main_out_dir')}\n\n"
    
    if 'receptor_pdb' in cfg or 'vina_receptor' in cfg or 'grid_receptor' in cfg:
        summary += f"{C.YELLOW}[Receptor & Search Space]{C.RESET}\n"
        if 'receptor_pdb' in cfg: summary += f"  Receptor PDB:  {cfg['receptor_pdb']}\n"
        if 'vina_receptor' in cfg: summary += f"  Receptor PDBQT:{cfg['vina_receptor']}\n"
        if 'grid_receptor' in cfg: summary += f"  Grid Receptor: {cfg['grid_receptor']}\n"
        
        if 'search_space' in cfg:
            summary += f"  Search Space:  {cfg['search_space'].capitalize()}\n"
            if cfg['search_space'] == 'targeted':
                if cfg.get('targeted_auto'):
                    summary += f"  Target Method: Autodetect from complex\n"
                else:
                    summary += f"  Target Method: Manual ({cfg.get('manual_res_clean')})\n"
        summary += "\n"
        
    if 'docking_type' in cfg:
        summary += f"{C.YELLOW}[Protein Preparation]{C.RESET}\n"
        summary += f"  Docking Type:  {cfg['docking_type'].capitalize()}\n"
        if cfg['docking_type'] == 'flexible':
            if cfg.get('flex_auto'):
                summary += f"  Flex Residues: Autodetect from complex ({cfg.get('prep_cocrystal', 'Same as Receptor')})\n"
            else:
                summary += f"  Flex Residues: Manual ({cfg.get('flex_res_clean')})\n"
        summary += "\n"
        
    if 'ligands_input' in cfg:
        summary += f"{C.YELLOW}[Ligand Preparation]{C.RESET}\n"
        summary += f"  Input Sources: {cfg['ligands_input']}\n\n"

    print(summary)
    choice = input(f"{C.GREEN}Are all these inputs correct? Proceed with run? (y/n): {C.RESET}").strip().lower()
    
    if choice != 'y':
        print(f"{C.RED}Run cancelled by user.{C.RESET}")
        sys.exit(0)
        
    os.makedirs(cfg['main_out_dir'], exist_ok=True)
    clean_summary = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])').sub('', summary)
    with open(cfg['log_file'], 'w', encoding='utf-8') as f:
        f.write("=================================================================\n")
        f.write("                   VINADOCKX INPUT LOG\n")
        f.write("=================================================================\n")
        f.write(f"Run ID: {cfg['run_id']}\n")
        f.write(clean_summary + "\n\n")

# ==========================================
# Pipeline Execution
# ==========================================
def run_pipeline(steps, cfg, paths):
    log_file = cfg['log_file']
    print_banner(log_file, to_terminal=False) 
    
    results = {"success": 0, "failed": 0, "logs": []}
    out_paths = {}

    def execute_script(script_key, args, step_desc, pbar, expected_out_dir=None, expected_ext=None):
        args = [str(a) if a is not None else "" for a in args]
        cmd = [sys.executable, "-u", paths[script_key]] + args
        
        log_msg(f"\n{C.BOLD}{C.MAGENTA}[>>>] Running: {script_key}{C.RESET}", log_file, use_tqdm=True)
        
        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"
        env["FORCE_COLOR"] = "1" 
        env["PYTHONUNBUFFERED"] = "1" 
        env["PYTHONUTF8"] = "1" 
        
        try:
            pbar.set_description(f"{C.MAGENTA}Running: {script_key:<13}{C.RESET}")
            
            process = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                env=env
            )
            
            sys.stdout.write(C.WHITE)
            sys.stdout.flush()
            
            while True:
                chunk = process.stdout.read1(1024) 
                if not chunk:
                    if process.poll() is not None:
                        break
                    continue
                
                text = chunk.decode('utf-8', errors='replace')
                sys.stdout.write(text)
                sys.stdout.flush()
                
                with open(log_file, 'a', encoding='utf-8') as f:
                    clean_text = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])').sub('', text)
                    f.write(clean_text)
                    
            sys.stdout.write(C.RESET)
            sys.stdout.flush()
            
            if process.returncode == 0:
                if expected_out_dir:
                    if not os.path.exists(expected_out_dir):
                        results['failed'] += 1
                        results['logs'].append(f"{script_key}: FAILED (Output directory missing)")
                        return False
                        
                    if expected_ext:
                        if isinstance(expected_ext, str):
                            expected_ext = (expected_ext,)
                        found_files = [f for f in os.listdir(expected_out_dir) if f.endswith(expected_ext)]
                    else:
                        # Allow finding ANY file/folder to account for VinaDock6's subdirectories
                        found_files = os.listdir(expected_out_dir)
                        
                    if not found_files:
                        results['failed'] += 1
                        ext_str = "/".join(expected_ext) if expected_ext else "outputs"
                        results['logs'].append(f"{script_key}: FAILED (No valid {ext_str} generated)")
                        return False

                results['success'] += 1
                results['logs'].append(f"{script_key}: SUCCESS - {step_desc}")
                return True
            else:
                results['failed'] += 1
                results['logs'].append(f"{script_key}: FAILED (Code {process.returncode})")
                return False
        except Exception as e:
            results['failed'] += 1
            results['logs'].append(f"{script_key}: ERROR - {str(e)}")
            return False

    print(f"\n{C.BOLD}{C.CYAN}Commencing Pipeline Execution...{C.RESET}\n")
    
    with tqdm(total=len(steps), bar_format="{l_bar}%s{bar}%s| {n_fmt}/{total_fmt}" % (Fore.LIGHTCYAN_EX, Fore.RESET), ncols=80) as pbar:

        # Step 1: InterResFi
        if 1 in steps:
            out_1 = os.path.join(cfg['main_out_dir'], "1 Interacting Residues - InterResFi")
            os.makedirs(out_1, exist_ok=True)
            pdb_in = cfg.get('cocrystal_pdb', cfg.get('receptor_pdb', ''))
            
            success = execute_script('InterResFi', ['-i', pdb_in, '-o', out_1, '-q'], "Interacting residues extracted", pbar, expected_out_dir=out_1, expected_ext=".txt")
            if success:
                txts = glob.glob(os.path.join(out_1, "*_InterRes.txt"))
                if txts: out_paths['inter_res_txt'] = txts[0]
            pbar.update(1)

        # Step 2: GridConfigGen
        if 2 in steps:
            cascade_fail = False
            out_2 = os.path.join(cfg['main_out_dir'], "2 Grid & Config - GridConfigGen")
            os.makedirs(out_2, exist_ok=True)
            args = ['--mode', cfg.get('search_space', 'blind'), '--pdb', cfg.get('receptor_pdb', ''), '--out', out_2, '--type', 'config']
            
            if cfg.get('search_space') == 'targeted':
                if cfg.get('targeted_auto'):
                    if 'inter_res_txt' in out_paths:
                        args.extend(['--residues', out_paths['inter_res_txt']])
                    else:
                        cascade_fail = True
                        log_msg(f"\n{C.RED}[!] Cascading Failure: No interacting residues found for GridConfigGen.{C.RESET}", log_file, use_tqdm=True)
                else:
                    args.extend(['--residues', cfg.get('manual_res_clean', '')])
                    
            if cascade_fail:
                results['failed'] += 1
                results['logs'].append("GridConfigGen: SKIPPED (Cascading Failure)")
                pbar.update(1)
            else:
                success = execute_script('GridConfigGen', args, "Grid Config Generated", pbar, expected_out_dir=out_2, expected_ext=".txt")
                if success:
                    txts = glob.glob(os.path.join(out_2, "*_Config_*.txt"))
                    if txts: out_paths['config_txt'] = txts[0]
                pbar.update(1)

        # Step 3: GridViz
        if 3 in steps:
            out_3 = os.path.join(cfg['main_out_dir'], "3 Grid Visualizer - GridViz")
            os.makedirs(out_3, exist_ok=True)
            cfg_txt = out_paths.get('config_txt', cfg.get('config_txt', ''))
            
            if not cfg_txt:
                log_msg(f"\n{C.RED}[!] Cascading Failure: No config file available for GridViz.{C.RESET}", log_file, use_tqdm=True)
                results['failed'] += 1
                results['logs'].append("GridViz: SKIPPED (Cascading Failure)")
                pbar.update(1)
            else:
                args = ['-c', cfg_txt, '-o', out_3]
                grid_rec = cfg.get('receptor_pdb') or cfg.get('vina_receptor') or cfg.get('grid_receptor')
                if grid_rec: args.extend(['-p', grid_rec])
                execute_script('GridViz', args, "Grid Visualized", pbar, expected_out_dir=out_3, expected_ext=None)
                pbar.update(1)

        # Step 4: PrepProt
        if 4 in steps:
            cascade_fail = False
            out_4 = os.path.join(cfg['main_out_dir'], "4 Protein Preparation - PrepProt")
            os.makedirs(out_4, exist_ok=True)
            args = ['-i', cfg.get('receptor_pdb', ''), '-o', out_4]
            
            if cfg.get('docking_type') == 'flexible':
                if cfg.get('flex_auto'):
                    if 'inter_res_txt' in out_paths:
                        flex_str = extract_residues(out_paths['inter_res_txt'])
                        if flex_str: 
                            args.extend(['-f', flex_str])
                        else:
                            cascade_fail = True
                    else:
                        out_tmp = os.path.join(cfg['main_out_dir'], "1 Interacting Residues - InterResFi")
                        os.makedirs(out_tmp, exist_ok=True)
                        execute_script('InterResFi', ['-i', cfg.get('prep_cocrystal', cfg.get('receptor_pdb', '')), '-o', out_tmp, '-q'], "Flex residues extracted", pbar, expected_out_dir=out_tmp, expected_ext=".txt")
                        txts = glob.glob(os.path.join(out_tmp, "*_InterRes.txt"))
                        if txts: 
                            flex_str = extract_residues(txts[0])
                            if flex_str: 
                                args.extend(['-f', flex_str])
                            else: cascade_fail = True
                        else: cascade_fail = True
                else:
                    args.extend(['-f', cfg.get('flex_res_clean', '')])
                    
            if cascade_fail:
                log_msg(f"\n{C.RED}[!] Cascading Failure: PrepProt failed to retrieve flexible residues.{C.RESET}", log_file, use_tqdm=True)
                results['failed'] += 1
                results['logs'].append("PrepProt: SKIPPED (Cascading Failure)")
                pbar.update(1)
            else:
                success = execute_script('PrepProt', args, "Protein Prepared", pbar, expected_out_dir=out_4, expected_ext=".pdbqt")
                if success: out_paths['prot_dir'] = out_4
                pbar.update(1)

        # Step 5: PrepLig
        if 5 in steps:
            out_5 = os.path.join(cfg['main_out_dir'], "5 Ligand Preparation - PrepLig")
            os.makedirs(out_5, exist_ok=True)
            success = execute_script('PrepLig', ['-i', cfg.get('ligands_input', ''), '-o', out_5, '-s', '1'], "Ligands Prepared", pbar, expected_out_dir=out_5, expected_ext=".pdbqt")
            if success: out_paths['lig_dir'] = out_5
            pbar.update(1)

        # Step 6: VinaDock
        if 6 in steps:
            cascade_fail = False
            out_6 = os.path.join(cfg['main_out_dir'], "6 Docking Results - VinaDock")
            tmp_dir = os.path.join(cfg['main_out_dir'], ".tmp_vina_inputs")
            os.makedirs(tmp_dir, exist_ok=True)
            
            rec_path, flex_path = "", ""
            
            if 4 in steps:
                pdbs = glob.glob(os.path.join(out_paths.get('prot_dir', ''), "*.pdbqt"))
                if not pdbs:
                    cascade_fail = True
                    log_msg(f"\n{C.RED}[!] Cascading Failure: PrepProt failed to generate a valid receptor PDBQT.{C.RESET}", log_file, use_tqdm=True)
                else:
                    for p in pdbs:
                        new_name = os.path.basename(p).replace('_PrepProt', '')
                        dest = os.path.join(tmp_dir, new_name)
                        shutil.copy(p, dest)
                        if '_flex' in new_name: flex_path = dest
                        elif '_rigid' in new_name or cfg.get('docking_type') == 'standard': rec_path = dest
            else:
                rec_path = cfg.get('vina_receptor', '')
                if not os.path.exists(rec_path):
                    cascade_fail = True
                    log_msg(f"\n{C.RED}[!] Cascading Failure: Provided receptor PDBQT file not found.{C.RESET}", log_file, use_tqdm=True)

            lig_src = out_paths.get('lig_dir', cfg.get('vina_ligands', ''))
            lig_tmp_dir = os.path.join(tmp_dir, "ligands")
            os.makedirs(lig_tmp_dir, exist_ok=True)
            valid_ligands_found = False
            
            if not cascade_fail:
                if os.path.isdir(lig_src):
                    for lf in glob.glob(os.path.join(lig_src, "*.pdbqt")):
                        shutil.copy(lf, os.path.join(lig_tmp_dir, os.path.basename(lf).replace('_PrepLig', '')))
                        valid_ligands_found = True
                else:
                    for lf in lig_src.split(','):
                        lf_clean = lf.strip()
                        if os.path.exists(lf_clean):
                            shutil.copy(lf_clean, os.path.join(lig_tmp_dir, os.path.basename(lf_clean).replace('_PrepLig', '')))
                            valid_ligands_found = True
                            
                if not valid_ligands_found:
                    cascade_fail = True
                    log_msg(f"\n{C.RED}[!] Cascading Failure: No valid ligand PDBQTs found to dock.{C.RESET}", log_file, use_tqdm=True)

            if cascade_fail:
                results['failed'] += 1
                results['logs'].append("VinaDock: SKIPPED (Cascading Failure - Missing inputs)")
                pbar.update(1)
            else:
                cfg_txt = out_paths.get('config_txt', cfg.get('config_txt', ''))
                args = ['--cli', '--outdir', out_6, '--config', cfg_txt, '--ligands', lig_tmp_dir]
                if flex_path: args.extend(['--receptor', rec_path, '--flex', flex_path])
                else: args.extend(['--receptor', rec_path])

                # Expected ext is None because VinaDock6.py puts things in subdirectories now
                execute_script('VinaDock', args, "Docking Completed", pbar, expected_out_dir=out_6, expected_ext=None)
                out_paths['vina_out'] = out_6
                pbar.update(1)

        # Step 7: BindReSort
        if 7 in steps:
            out_7 = os.path.join(cfg['main_out_dir'], "7 Sorted Results - BindReSort")
            in_7 = out_paths.get('vina_out', cfg.get('bindresort_input', ''))
            
            valid_logs = False
            if in_7 and os.path.exists(in_7):
                # Search recursively to find logs inside the new "Log" subdirectory created by VinaDock6
                if os.path.isdir(in_7) and glob.glob(os.path.join(in_7, "**", "*.txt"), recursive=True): 
                    valid_logs = True
                elif os.path.isfile(in_7) and in_7.endswith('.txt'): 
                    valid_logs = True
                
            if not valid_logs:
                log_msg(f"\n{C.RED}[!] Cascading Failure: No Vina log files found for BindReSort.{C.RESET}", log_file, use_tqdm=True)
                results['failed'] += 1
                results['logs'].append("BindReSort: SKIPPED (Cascading Failure - No logs to sort)")
                pbar.update(1)
            else:
                execute_script('BindReSort', ['-i', in_7, '-o', out_7, '-f', '2'], "Results Sorted", pbar, expected_out_dir=out_7, expected_ext=".xlsx")
                pbar.update(1)
                
        pbar.set_description(f"{C.GREEN}Pipeline Complete!{C.RESET}")

    # --- MASTER CLEANUP ---
    tmp_dir_path = os.path.join(cfg['main_out_dir'], ".tmp_vina_inputs")
    if os.path.exists(tmp_dir_path):
        time.sleep(1) 
        shutil.rmtree(tmp_dir_path, ignore_errors=True)

    # Final Detailed Summary
    summary = f"""
{C.BOLD}{C.MAGENTA}================================================================={C.RESET}
{C.BOLD}{C.WHITE}                          FINAL RUN SUMMARY{C.RESET}
{C.BOLD}{C.MAGENTA}================================================================={C.RESET}
{C.WHITE}Total Scripts Run: {C.CYAN}{len(steps)}{C.RESET}
{C.WHITE}Successful Runs:   {C.GREEN}{results['success']}{C.RESET}
{C.WHITE}Failed Runs:       {C.RED}{results['failed']}{C.RESET}

{C.CYAN}Master Output Directory:{C.RESET} {os.path.abspath(cfg['main_out_dir'])}
{C.CYAN}Master Log File Saved:  {C.RESET} {os.path.basename(cfg['log_file'])}

{C.YELLOW}Detailed Sub-script Executions:{C.RESET}
"""
    for log in results['logs']:
        if "SUCCESS" in log: summary += f"  {C.GREEN}[✅]{C.RESET} {log}\n"
        elif "SKIPPED" in log: summary += f"  {C.YELLOW}[⭕]{C.RESET} {log}\n"
        else: summary += f"  {C.RED}[❌]{C.RESET} {log}\n"
            
    summary += f"{C.MAGENTA}================================================================={C.RESET}\n"
    log_msg(summary, log_file, use_tqdm=True)

# ==========================================
# Main Menu
# ==========================================
def main():
    paths = find_sub_scripts()
    print_banner()
    
    print(f"{C.BOLD}{C.MAGENTA}--- MAIN MENU ---{C.RESET}")
    print(f"  {C.CYAN}1.{C.RESET} Complete Docking Pipeline (Scripts 1-7)")
    print(f"  {C.CYAN}2.{C.RESET} Run Pipeline Partially")
    choice = clean_input("Select an option (1 or 2):")
    
    if choice == '1':
        steps = list(range(1, 8))
    else:
        print(f"\n{C.BOLD}Available Pipeline Scripts:{C.RESET}")
        print(f"  {C.CYAN}1{C.RESET} InterResFi: Find Interacting Residues")
        print(f"  {C.CYAN}2{C.RESET} GridConfigGen: Generate Grid & Config")
        print(f"  {C.CYAN}3{C.RESET} GridViz: Visualize grid")
        print(f"  {C.CYAN}4{C.RESET} PrepProt: Protein preparation")
        print(f"  {C.CYAN}5{C.RESET} PrepLig: Ligand Preparation")
        print(f"  {C.CYAN}6{C.RESET} VinaDock: Docking Execution")
        print(f"  {C.CYAN}7{C.RESET} BindReSort: Sort docking results")
        steps = parse_run_sequence(clean_input("Enter sequence:", "1, 3, 5-7"))
        
    check_pipeline_dependencies(steps)
    cfg = gather_inputs(steps)
    confirm_inputs(cfg, steps)
    run_pipeline(steps, cfg, paths)

if __name__ == "__main__":
    main()