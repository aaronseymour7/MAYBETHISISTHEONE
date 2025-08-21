#!/usr/bin/env python3

import glob
import os
import csv
from hfrpkg.compute_folder import compute_folder
from hfrpkg.write_data import write_single_reaction

def main():
    results = []

    with open("Reactions_summary.txt", "w", encoding="utf-8") as rxn_fout:
        for mhfr_file in sorted(glob.glob("*.mhfr"), key=lambda x: int(x.split(".")[0])):
            if os.path.isdir(mhfr_file):
                try:
                    data = compute_folder(mhfr_file)
                    if data is not None:
                        result = write_single_reaction(data, mhfr_file)
                        if result is not None:
                            results.append(result)

                            summary_file = os.path.join(mhfr_file, "opt_summary.txt")
                            with open(summary_file, "r", encoding="utf-8") as sf:
                                rxn_fout.write(f"=== {mhfr_file}/reaction_summary.txt ===\n")
                                rxn_fout.write(sf.read())
                                rxn_fout.write("\n")

                except Exception as e:
                    print(f"Error processing {mhfr_file}: {e}")

    with open("enthalpies_summary.csv", "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(["SMILES", "LEVEL", "InChI", "ΔHf DFT (kcal/mol)", "ΔHf ATcT (kcal/mol)"])
        writer.writerows(results)

if __name__ == "__main__":
    main()

