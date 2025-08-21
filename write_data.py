import os
import glob
import re
from hfrpkg.compute_folder import compute_folder

def kjTokcal(value):
        try:
            return value / 4.184
        except TypeError:
            return None  
def write_single_reaction(reaction_data, folder_path):
    output_path = os.path.join(folder_path, "opt_summary.txt")
    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write(f"{reaction_data['input_smiles']}\t{reaction_data['input_inchi']}\n")
        fout.write(f"Reaction Enthalpy (kJ/mol):{reaction_data['reaction_H']}\n")
        fout.write(f"DFT Enthalpy of Formation (kJ/mol):{reaction_data['dft_hf']}\n")
        fout.write("REACTANTS\n")
        for coeff, smiles, inchi, atct, energy, zpve in reaction_data['reactants']:
            fout.write(f"{coeff} {smiles}\t{inchi}\t{atct}\t{energy}\t{zpve}\n")
        
        fout.write("PRODUCTS\n")
        for coeff, smiles, inchi, atct, energy, zpve in reaction_data['products']:
            fout.write(f"{coeff} {smiles}\t{inchi}\t{atct}\t{energy}\t{zpve}\n")
        
        fout.write("\n")
    

    return (
        reaction_data['input_smiles'],
        reaction_data['level'],
        reaction_data['input_inchi'],
        kjTokcal(reaction_data['dft_hf']),
        kjTokcal(reaction_data['input_atct']),
    )
