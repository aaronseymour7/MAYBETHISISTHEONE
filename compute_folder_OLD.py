#!/usr/bin/env python3

import os
import glob
import re
import sys
from AaronTools.fileIO import FileReader

def run_folder(folder_path):
    import os
    import glob
    import re
    from AaronTools.fileIO import FileReader

    def get_B(filename):
        base = filename.replace(".log", "")
        return int(base.split("_")[0][1:])

    def get_final_R_coefficient():
        pattern = re.compile(r"^R(\d+)_(\d+)\.log$")
        max_B = -1
        coeff_C = None
        for filename in glob.glob("*.log"):
            match = pattern.match(filename)
            if match:
                B, C = int(match.group(1)), int(match.group(2))
                if B > max_B:
                    max_B, coeff_C = B, C
        return coeff_C

    def extract_coeff_and_type(filename):
        m = re.match(r"^([PR])\d+_(\d+)\.log$", filename)
        if not m:
            return None, None
        return filename.replace(".log", ""), int(m.group(2))

    def get_enthalpy(logfile):
        try:
            reader = FileReader(logfile, just_geom=False)
            if 'enthalpy' in reader.keys():
                return reader['enthalpy']
            else:
                return None
        except Exception:
            return None

    def get_inchi(log_filename, index_path="index.txt"):
        try:
            with open(index_path) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2 and parts[0] == log_filename:
                        return parts[1]
        except FileNotFoundError:
            pass
        return None

    def get_Hf(inchi):
        #atct_path="/home/ads09449/bin/ATcT_lib.txt"
        try:
            with importlib.resources.open_text("hfrpkg.data", "ATcT_lib.txt", encoding="utf-8") as f:
            #with open(atct_path, encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6 and parts[3] == inchi:
                        return float(parts[5])
        except Exception:
            pass
        return None

    cwd = os.getcwd()
    os.chdir(folder_path)

    try:
        log_files = sorted(glob.glob("*.log"), key=get_B)
        products = [f for f in log_files if f.startswith("P")]
        reactants = [f for f in log_files if f.startswith("R")]

        if not reactants:
            raise ValueError("No reactants found.")

        total_products = 0.0
        total_reactants = 0.0
        Hf_products = 0.0
        Hf_reactants = 0.0
        missing_Hf = False

        for f in products:
            mol_type, coeff = extract_coeff_and_type(f)
            if mol_type is None: continue
            enthalpy = get_enthalpy(f)
            if enthalpy is not None:
                total_products += coeff * enthalpy

        for f in reactants:
            mol_type, coeff = extract_coeff_and_type(f)
            if mol_type is None: continue
            enthalpy = get_enthalpy(f)
            if enthalpy is not None:
                total_reactants += coeff * enthalpy

        reaction = 2625.5 * (total_products - total_reactants)
        final_coeff = get_final_R_coefficient()
        if final_coeff is None:
            raise ValueError("Could not determine final coefficient")

        for f in reactants[:-1]:  # usable reactants
            mol_type, coeff = extract_coeff_and_type(f)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            Hf = get_Hf(inchi)
            if Hf is None:
                print(f"[ATcT MISSING]  {folder_path}  {mol_type} → InChI: {inchi}")
                missing_Hf = True
            else:
                Hf_reactants += Hf * coeff

        for f in products:
            mol_type, coeff = extract_coeff_and_type(f)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            Hf = get_Hf(inchi)
            if Hf is None:
                print(f"[ATcT MISSING]  {folder_path} {mol_type} → InChI: {inchi}")
                missing_Hf = True
            else:
                Hf_products += Hf * coeff

        if missing_Hf:
            input_hf = ""
        else:
            input_hf = round((Hf_products - Hf_reactants - reaction) / (4.184*final_coeff), 2)

        target_file = reactants[-1]
        mol_type, _ = extract_coeff_and_type(target_file)
        inchi = get_inchi(mol_type)
        if get_Hf(inchi):
            atct = get_Hf(inchi)/4.184
        else:
            atct = get_Hf(inchi)
        atct_str = round(atct, 2) if atct is not None else ""
        smiles = ""
        try:
            with open("index.txt") as f:
                lines = f.readlines()
                if lines:
                    header = lines[0].strip().split("\t")
                    level = header[1]
                for line in lines[2:]:
                    parts = line.strip().split("\t")
                    if parts[0] == mol_type and len(parts) >= 3:
                        smiles = parts[2]
                        break
        except FileNotFoundError:
            print(f"index.txt not found in {folder_path}")
        return smiles, level, inchi, input_hf, atct_str
    except Exception as e:
        print(f"[ERROR] Folder {folder_path}: {e}")
        return folder_path, "ERROR", ""

    finally:
        os.chdir(cwd)
