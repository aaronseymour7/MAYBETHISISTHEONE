    #!/usr/bin/env python3

import os
import glob
import re
from AaronTools.fileIO import FileReader
import importlib.resources  

def compute_folder(folder_path):
    def get_B(filename):
        base = os.path.splitext(filename)[0]
        return int(base.split("_")[0][1:])

    def get_final_R_coefficient(out=".log"):
        pattern = re.compile(rf"^R(\d+)_(\d+){re.escape(out)}$")
        max_B = -1
        coeff_C = None
        for filename in glob.glob("*" + out):
            match = pattern.match(filename)
            if match:
                B, C = int(match.group(1)), int(match.group(2))
                if B > max_B:
                    max_B, coeff_C = B, C
        return coeff_C

    def extract_coeff_and_type(filename, out=".log"):
        escaped_ext = re.escape(out)  # in case the dot needs escaping
        m = re.match(rf"^([PR])\d+_(\d+){escaped_ext}$", filename)
        if not m:
            return None, None
        return os.path.splitext(filename)[0], int(m.group(2))

    def get_enthalpy(logfile):
        try:
            reader = FileReader(logfile, just_geom=False)
            if 'E_ZPVE' in reader.keys():
                return reader['E_ZPVE']
            else:
                return None
        except Exception:
            return None
    def get_zpve(logfile):
        try:
            reader = FileReader(logfile, just_geom=False)
            if 'ZPVE' in reader.keys():
                return reader['ZPVE']
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

    def get_smiles(log_filename, index_path="index.txt"):
        try:
            with open(index_path) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 3 and parts[0] == log_filename:
                        return parts[2]
        except FileNotFoundError:
            pass
        return None

    def get_level(index_path="index.txt"):
        try:
            with open(index_path) as f:
                first_line = f.readline()
                # Expecting format: "Level:\t reaction_fn_name"
                parts = first_line.strip().split("\t")
                if len(parts) >= 2:
                    return parts[1]
        except FileNotFoundError:
            pass
        return None

    def get_Hf(inchi):
        try:
            with importlib.resources.open_text("hfrpkg.data", "ATcT_lib.txt", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6 and parts[3] == inchi:
                        return float(parts[5])
        except Exception:
            pass
        return None
    
    cwd = os.getcwd()
    os.chdir(folder_path)
    
    def get_extension(index_path="index.txt"):
        ext_map = {
            "gaussian": ".log",
            "orca": ".out",
            "psi4": ".dat"
        }
        try:
            with open(index_path) as f:
                first_line = f.readline()
                # Expecting format: "Level:\t reaction_fn_name\t software: \t software_name"
                parts = first_line.strip().split("\t")
                if len(parts) >= 4:
                    ext = ext_map.get(parts[3].lower())
                    if ext is None:
                        print(f"[ERROR] Unknown software '{parts[3]}' in index.txt.")
                    return ext
        except FileNotFoundError:
            pass
        return None
    
    ext = get_extension()


    try:
        log_files = sorted(glob.glob("*"+ ext), key=get_B)
        products = [f for f in log_files if f.startswith("P")]
        reactants = [f for f in log_files if f.startswith("R")]

        if not reactants:
            raise ValueError("No reactants found.")

        total_products = 0.0
        total_reactants = 0.0
        Hf_products = 0.0
        Hf_reactants = 0.0
        missing_Hf = False

        reactants_data = []
        products_data = []

        # Sum DFT enthalpies for reactants/products
        for f in products:
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            enthalpy = get_enthalpy(f)
            if enthalpy is not None:
                total_products += coeff * enthalpy

        for f in reactants:
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            enthalpy = get_enthalpy(f)
            if enthalpy is not None:
                total_reactants += coeff * enthalpy

        reaction = 2625.5 * (total_products - total_reactants)
        final_coeff = get_final_R_coefficient(ext)
        if final_coeff is None:
            raise ValueError("Could not determine final coefficient")

        # Collect reactant info
        for f in reactants[:-1]:  
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            smiles = get_smiles(mol_type)
            Hf = get_Hf(inchi)
            energy = get_enthalpy(f)
            zpve = get_zpve(f)
            if Hf is None:
                print(f"[ATcT MISSING]  {folder_path}  {mol_type} → InChI: {inchi}")
                missing_Hf = True
                Hf = ""
            Hf_reactants += (Hf if Hf else 0) * coeff
            reactants_data.append((coeff, smiles, inchi, Hf, energy, zpve))
        input_file = reactants[-1]
        mol_type, coeff = extract_coeff_and_type(input_file, ext)
        if mol_type is not None:
            inchi = get_inchi(mol_type)
            input_inchi = inchi
            smiles = get_smiles(mol_type)
            input_smiles = smiles
            Hf = get_Hf(inchi)
            energy = get_enthalpy(input_file)
            zpve = get_zpve(input_file)
            if Hf is None:
                print(f"[ATcT MISSING]  {folder_path}  {mol_type} → InChI: {inchi}")
                Hf = ""
        atct_value = Hf
        reactants_data.append((coeff, smiles, inchi, Hf, energy, zpve))
        
        # Collect product info
        for f in products:
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            smiles = get_smiles(mol_type)
            Hf = get_Hf(inchi)
            energy = get_enthalpy(f)
            zpve = get_zpve(f)
            if Hf is None:
                print(f"[ATcT MISSING]  {folder_path} {mol_type} → InChI: {inchi}")
                missing_Hf = True
                Hf = ""
            Hf_products += (Hf if Hf else 0) * coeff
            products_data.append((coeff, smiles, inchi, Hf, energy, zpve))


        if missing_Hf:
            input_hf = ""
        else:
            input_hf = round((Hf_products - Hf_reactants - reaction) / (final_coeff), 4)
    
        

        level = get_level()

        return {
            "input_smiles": input_smiles,
            "input_inchi": input_inchi,
            "dft_hf": input_hf,
            "level": level,
            "reactants": reactants_data,
            "products": products_data,
            "input_atct": atct_value,
            "reaction_H": reaction
        }

    except Exception as e:
        print(f"[ERROR] Folder {folder_path}: {e}")
        return None

    finally:
        os.chdir(cwd)
