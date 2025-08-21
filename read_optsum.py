import os

def read_optsum(folder_path):

    input_path = os.path.join(folder_path, "opt_summary.txt")
    reaction_data = {
        "input_smiles": None,
        "input_inchi": None,
        "reaction_H": None,
        "dft_hf": None,
        "reactants": [],
        "products": [],
    }

    with open(input_path, "r", encoding="utf-8") as fin:
        lines = [line.strip() for line in fin if line.strip()]

    first_parts = lines[0].split("\t")
    reaction_data["input_smiles"] = first_parts[0]
    reaction_data["input_inchi"] = first_parts[1] if len(first_parts) > 1 else None

    reaction_data["reaction_H"] = float(lines[1].split(":")[1])
   # reaction_data["dft_hf"] = float(lines[2].split(":")[1])
    
    try:
        value = lines[2].split(":")[1].strip()
        reaction_data["dft_hf"] = float(value) if value else None
    except (IndexError, ValueError):
        reaction_data["dft_hf"] = None
    
    section = None
    for line in lines[3:]:
        if line == "REACTANTS":
            section = "reactants"
            continue
        elif line == "PRODUCTS":
            section = "products"
            continue

        if section:
            parts = line.split("\t")
            coeff_smiles = parts[0].split(" ")
            coeff = int(coeff_smiles[0])
            smiles = coeff_smiles[1]
            inchi = parts[1]
           # atct = float(parts[2]) if parts[2] != "None" else None
            try:
                value = parts[2]
                atct = float(value) if value else None
            except (IndexError, ValueError):
                atct = None
            energy = float(parts[3]) if parts[3] != "None" else None
            zpve = float(parts[4]) if parts[4] != "None" else None

            reaction_data[section].append((coeff, smiles, inchi, atct, energy, zpve))

    return reaction_data


