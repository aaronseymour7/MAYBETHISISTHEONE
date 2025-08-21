    #!/usr/bin/env python3

import os
import glob
import argparse
import re
from AaronTools.fileIO import FileReader
import importlib.resources  
from hfrpkg.read_optsum import read_optsum
from hfrpkg.compute_folder import compute_folder
from hfrpkg.write_data import write_single_reaction


def spec_compute(mhfr_file=None):
    notSP = False
    cwd = os.getcwd()  
    if mhfr_file is None:
        notSP = True
        mhfr_file = cwd
    
    reaction_data = compute_folder(mhfr_file)
    write_single_reaction(reaction_data, mhfr_file)  
    opt_data = read_optsum(mhfr_file)                

    
    spec_dir = os.path.join(mhfr_file, "spec")
    os.chdir(spec_dir)
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

    def get_enthalpy(logfile, inchi):
        zpve = None
        for coeff, smiles, r_inchi, atct, energy, zpve_val in (opt_data["reactants"] + opt_data["products"]):
            if r_inchi == inchi:
                zpve = zpve_val
                break
        

        try:
            reader = FileReader(logfile, just_geom=False)
            if "energy" in reader.keys():
                if zpve is None:
                    print(f"No ZPVE found in opt_summary for {inchi}")
                    return None, zpve
                return reader["energy"] + zpve, zpve
        except Exception:
            pass
        return None, None
    

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
    
    
    
    def get_extension(index_path="index.txt"):
        ext_map = {
            "gaussian": ".log",
            "orca": ".out",
            "psi4": ".out"
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
    #print(f"Using {ext} files...")

    try:
        log_files = sorted(glob.glob("*"+ ext), key=get_B)
        #if log_files:
            #print("spoutputs found")
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

        for f in products:
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            enthalpy, zpve = get_enthalpy(f, inchi)
            if enthalpy is not None:
                total_products += coeff * enthalpy

        for f in reactants:
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            enthalpy, zpve = get_enthalpy(f, inchi)
            if enthalpy is not None:
                total_reactants += coeff * enthalpy

        reaction = 2625.5 * (total_products - total_reactants)
        final_coeff = get_final_R_coefficient(ext)
        if final_coeff is None:
            raise ValueError("Could not determine final coefficient")

        for f in reactants[:-1]:  
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            smiles = get_smiles(f)
            Hf = get_Hf(inchi)
            enthalpy, zpve = get_enthalpy(f, inchi)
            
            if Hf is None:
                print(f"[ATcT MISSING]  {spec_dir}  {mol_type} → InChI: {inchi}")
                missing_Hf = True
            Hf_reactants += (Hf if Hf else 0) * coeff
            reactants_data.append((coeff, smiles, inchi, Hf, enthalpy, zpve))
        input_file = reactants[-1]
        mol_type, coeff = extract_coeff_and_type(input_file, ext)
        if mol_type is not None:
            inchi = get_inchi(mol_type)
            input_smiles = get_smiles(input_file)
            input_inchi = inchi        
            Hf = get_Hf(inchi)
            enthalpy, zpve = get_enthalpy(input_file, inchi)
            
            if Hf is None:
                print(f"[ATcT MISSING]  {spec_dir}  {mol_type} → InChI: {inchi}")
            atct_value = Hf
            reactants_data.append((coeff, input_smiles, inchi, Hf, enthalpy, zpve))
        
        for f in products:
            mol_type, coeff = extract_coeff_and_type(f, ext)
            if mol_type is None: continue
            inchi = get_inchi(mol_type)
            smiles = get_smiles(f)
            Hf = get_Hf(inchi)
            enthalpy, zpve = get_enthalpy(f, inchi)
            if Hf is None:
                print(f"[ATcT MISSING]  {spec_dir} {mol_type} → InChI: {inchi}")
                missing_Hf = True
            Hf_products += (Hf if Hf else 0) * coeff
            products_data.append((coeff, smiles, inchi, Hf, enthalpy, zpve))


        if missing_Hf:
            input_hf = None
        else:
            input_hf = round((Hf_products - Hf_reactants - reaction) / (final_coeff), 4)
    
        

        level = get_level()

        sp_data = {
            "input_smiles": input_smiles,
            "input_inchi": input_inchi,
            "dft_hf": input_hf,
            "level": level,
            "reactants": reactants_data,
            "products": products_data,
            "input_atct": atct_value,
            "reaction_H": reaction
        }
        #print(spec_dir)
        
        output_path = os.path.join(spec_dir, "sp_summary.txt")
        #print(output_path)
        with open(output_path, "w", encoding="utf-8") as fout:
            #print(output_path)
            fout.write(f"{sp_data['input_smiles']}\t{sp_data['input_inchi']}\n")
            fout.write(f"Reaction Enthalpy (kJ/mol):{sp_data['reaction_H']}\n")
            fout.write(f"DFT Enthalpy of Formation (kJ/mol):{sp_data['dft_hf']}\n")
            fout.write("REACTANTS\n")
            for coeff, smiles, inchi, atct, enthalpy, zpve in sp_data['reactants']:
                fout.write(f"{coeff} {smiles}\t{inchi}\t{atct}\t{enthalpy}\t{zpve}\n")
        
            fout.write("PRODUCTS\n")
            for coeff, smiles, inchi, atct, enthalpy, zpve in sp_data['products']:
                fout.write(f"{coeff} {smiles}\t{inchi}\t{atct}\t{enthalpy}\t{zpve}\n")
        #print("yeah this should work")
        return sp_data

    except Exception as e:
        print(f"[ERROR] Folder {spec_dir}: {e}")
        return None

    finally:
        os.chdir(cwd)

def main_cli():
    parser = argparse.ArgumentParser(description="Generate and submit SP .com files from optimized .log files")
    parser.add_argument("--f", "--folder", dest="f", default= None, help="mhfr folder, if only single, leave empty")
    args = parser.parse_args()

    spec_compute(mhfr_file=args.f)

if __name__ == "__main__":
   main_cli()
