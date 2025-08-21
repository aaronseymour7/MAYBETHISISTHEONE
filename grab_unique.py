
import os
import glob
import shutil
from hfrpkg.utils import get_extensions

def grab_unique_coms(software):
    output_dir = "unique_files"
    os.makedirs(output_dir, exist_ok=True)

    seen_inchis = {}
    counter = 1  

    output_index_path = os.path.join(output_dir, "index.txt")
    with open(output_index_path, "w") as index_out:
        index_out.write(f"Software:\t{software}\n")
        index_out.write("Filename\tInChI\tSMILES\n")
        for mhfr_dir in sorted(glob.glob("*.mhfr")):
            index_file_path = os.path.join(mhfr_dir, "index.txt")
            
            if not os.path.isfile(index_file_path):
                continue
            inext, outext = get_extensions(index_file_path)
            with open(index_file_path, "r") as f:
                lines = f.readlines()

            for line in lines:
                
                if line.startswith("Filename") or line.startswith("Level:"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue

                filename, inchi, smiles = parts
                if inchi in seen_inchis:
                    continue  
                
                source_com_path = os.path.join(mhfr_dir, filename + inext)
                if not os.path.isfile(source_com_path):
                    print(f"[WARNING] Missing .com file: {source_com_path}")
                    continue

                dest_com_filename = f"{counter}{inext}"
                dest_com_path = os.path.join(output_dir, dest_com_filename)

                shutil.copyfile(source_com_path, dest_com_path)

                index_out.write(f"{counter}\t{inchi}\t{smiles}\n")

                seen_inchis[inchi] = counter
                counter += 1

    print(f"Done. Unique .com files written to: {output_dir}/")
