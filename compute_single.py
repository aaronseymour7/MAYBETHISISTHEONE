#!/usr/bin/env python3

import glob
import os
import csv
from hfrpkg.compute_folder import compute_folder
from hfrpkg.write_data import write_single_reaction

def main():
    result = None
    mhfr_file = "."  

    try:
        data = compute_folder(mhfr_file)
        if data is not None:
            result = write_single_reaction(data, mhfr_file)
    except Exception as e:
        print(f"Error processing {mhfr_file}: {e}")

    with open("enthalpy_summary.csv", "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(["SMILES", "LEVEL", "InChI", "ΔHf DFT (kcal/mol)", "ΔHf ATcT (kcal/mol)"])
        writer.writerow(result)

def main_cli():
    main()


if __name__ == "__main__":
    main_cli()
