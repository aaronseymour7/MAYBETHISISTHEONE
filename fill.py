#!/usr/bin/env python3

import os
import glob
import shutil
from hfrpkg.utils import get_extensions
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
    unique_index_path = os.path.join(unique_folder, "index.txt")
    if not os.path.exists(unique_index_path):
        print("Missing unique_files/index.txt")
        return

    inchi_map = load_unique_inchi_map(unique_index_path)

    for mhfr_dir in glob.glob("*.mhfr"):
        index_file = os.path.join(mhfr_dir, "index.txt")
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
            src_log_path = os.path.join(unique_folder, unique_name)
            dst_log_path = os.path.join(mhfr_dir, mhfr_filename + outext)

            if not os.path.exists(src_log_path):
                print(f"[WARNING] Missing log file: {src_log_path}")
                continue

            shutil.copyfile(src_log_path, dst_log_path)
            #print(f"Copied {src_log_path} → {dst_log_path}")



