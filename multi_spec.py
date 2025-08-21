import os
import glob
import shutil
import csv
from hfrpkg.utils import get_extensions, get_softext_UFI
from hfrpkg.spec_compute import spec_compute


def load_unique_inchi_map(index_path):
    inchi_to_filename = {}
    with open(index_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2 and not parts[0].startswith("Level"):
                filename, inchi = parts[0], parts[1]
                inchi_to_filename[inchi] = filename
    return inchi_to_filename

def get_ext_from_soft(software):
    ext_map = {
        "gaussian": ('.com','.log'),
        "orca": ('.inp', '.out'),
        "psi4": ('.in', '.dat')
    }
    return ext_map.get(software.lower())

def fill_logs(unique_folder):
    unique_index_path = os.path.join(unique_folder, "spec", "index.txt")
    if not os.path.exists(unique_index_path):
        print("Missing unique_files/spec/index.txt")
        return
    software, x, y = get_softext_UFI(unique_index_path)
    inchi_map = load_unique_inchi_map(unique_index_path)

    for mhfr_dir in glob.glob("*.mhfr"):
        spec_dir = os.path.join(mhfr_dir, "spec")
        os.makedirs(spec_dir, exist_ok=True)

        index_path = os.path.join(mhfr_dir, "index.txt")
        if os.path.exists(index_path):
            with open(index_path, "r", encoding="utf-8") as f:
                lines = f.readlines()
            if lines and lines[0].startswith("Level:"):
                parts = lines[0].split("\t")
                if "Software:" in parts:
                    sw_idx = parts.index("Software:") + 1
                    if sw_idx < len(parts):
                        parts[sw_idx] = software
                        lines[0] = "\t".join(parts)+ "\n"

        with open(os.path.join(spec_dir, "index.txt"), "w", encoding="utf-8") as f:
            f.writelines(lines)

        index_file = os.path.join(mhfr_dir, "spec", "index.txt")
        if not os.path.exists(index_file):
            print(f"Skipping {mhfr_dir}, missing index.txt")
            continue

        inext, outext = get_extensions(index_file)
        with open(index_file, "r") as f:
            lines = f.readlines()[3:]

        for line in lines:
            if line.startswith("Level") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            mhfr_filename, inchi = parts[0], parts[1]
            if inchi not in inchi_map:
                print(f"[WARNING] InChI not found: {inchi}")
                continue

            unique_name = inchi_map[inchi] + outext
            src_log_path = os.path.join(unique_folder, "spec", unique_name)
            dst_log_path = os.path.join(mhfr_dir, "spec", mhfr_filename + outext)

            if not os.path.exists(src_log_path):
                print(f"[WARNING] Missing log file: {src_log_path}")
                continue

            os.makedirs(os.path.join(mhfr_dir, "spec"), exist_ok=True)
            shutil.copyfile(src_log_path, dst_log_path)

def main():
    fill_logs("unique_files")
    cwd = os.getcwd()
    results = []
    mhfrs = []
    for f in os.listdir():
        if f.endswith(".mhfr") and os.path.isdir(f):
            mhfrs.append(f)


    for mhfr_dir in mhfrs:
        try:
            os.chdir(mhfr_dir)
            sp_data = spec_compute()  

            if sp_data is not None:
                #print("greatnews")
                #results.append(sp_data)
                results.append([
                    sp_data["input_inchi"],
                    sp_data["level"],
                    sp_data["dft_hf"],      
                    sp_data["input_atct"]   
                ])
                summary_file = os.path.join(mhfr_dir, "spec", "sp_summary.txt")
                if os.path.exists(summary_file):
                    with open(summary_file, "r", encoding="utf-8") as sf:
                        rxn_fout.write(f"=== {mhfr_dir}/reaction_summary.txt ===\n")
                        rxn_fout.write(sf.read())
                        rxn_fout.write("\n")
            os.chdir(cwd)
        except Exception as e:
            print(f"Error processing {mhfr_dir}: {e}")
            os.chdir(cwd)
    '''
    with open("sp_reaction_summaries.txt", "w", encoding="utf-8") as rxn_fout:
        
        for mhfr_dir in logs:
            try:
                #os.chdir(mhfr_dir)
                sp_data = spec_compute(mhfr_dir)  # make sure this exists
                #os.chdir(cwd)

                if sp_data is not None:
                    results.append(sp_data)
                    summary_file = os.path.join(mhfr_dir, "sp_summary.txt")
                    if os.path.exists(summary_file):
                        with open(summary_file, "r", encoding="utf-8") as sf:
                            rxn_fout.write(f"=== {mhfr_dir}/reaction_summary.txt ===\n")
                            rxn_fout.write(sf.read())
                            rxn_fout.write("\n")
            except Exception as e:
                print(f"Error processing {mhfr_dir}: {e}")
                os.chdir(cwd)
    '''
    # write enthalpies summary
    with open("enthalpies_summary.csv", "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(["InChI", "LEVEL", "ΔHf DFT (kcal/mol)", "ΔHf ATcT (kcal/mol)"])
        writer.writerows(results)
    print("Summaries written to 'enthalpies_summary.csv'")
def main_cli():
    main()


if __name__ == "__main__":
    main_cli()
