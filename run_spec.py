

import sys
import os
import glob
import argparse
import shutil
from AaronTools.geometry import Geometry
from AaronTools.fileIO import FileWriter
from AaronTools.theory import Theory
from AaronTools.job_control import SubmitProcess
from AaronTools.theory.job_types import SinglePointJob
from hfrpkg.run_single import run_jobs
from hfrpkg.utils import get_extensions
def make_spec(method, basis, extension):
    
    
    folder = os.getcwd()
    optin, optout = get_extensions()
    log_files = glob.glob(os.path.join(folder, "*"+ optout))
    extension_map = {
        "g": ".com",
        "o": ".inp",
        "p": ".in"
    }
    software_map = {
        "g": "Gaussian",
        "o": "ORCA",
        "p": "Psi4"
    }
    spext = extension_map.get(extension.lower())
    software = software_map.get(extension.lower())
    if spext is None:
        print(f"Unknown software code '{extension}'. Please use 'g', 'o', or 'p'.")
        sys.exit(1)
    if not log_files:
        print("No .log files found in folder.")
        
        sys.exit(1)

    level = Theory(
        method=method,
        basis=basis,
        job_type=SinglePointJob()
    )
    com_files = []
    spec_dir = os.path.join(folder, "spec")
    os.makedirs(spec_dir, exist_ok=True)

    index_path = os.path.join(folder, "index.txt")
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

        print(f"Copied and updated index.txt to {spec_dir}")
    else:
        print("Warning: index.txt not found in parent directory.")
        
    for log_path in log_files:
        name = os.path.splitext(os.path.basename(log_path))[0]
        com_path = os.path.join(spec_dir, name + spext)

        try:
            geom = Geometry(log_path)
            geom.write(outfile=com_path, theory=level)
            print(f"Wrote: {com_path}")
        except Exception as e:
            print(f"Error processing {log_path}: {e}")
            
    cwd = os.getcwd()
    os.chdir(spec_dir)
    try:
        run_jobs()
    finally:
        os.chdir(cwd)

    
def main_cli():
    parser = argparse.ArgumentParser(description="Generate and submit SP .com files from optimized .log files")
    parser.add_argument("--m", "--method", dest="method", default="m06-2x", help="DFT method (default: m06-2x)")
    parser.add_argument("--b", "--basis", dest="basis", default="def2tzvp", help="Basis set (default: def2tzvp)")
    parser.add_argument(
        "--s", choices=["g", "o", "p"], default="g",
        help="Software to use: 'g' for Gaussian, 'o' for Orca, 'p' for Psi4 (default: 'g')"
    )
    args = parser.parse_args()

    make_spec(method=args.method, basis=args.basis, extension = args.s)

if __name__ == "__main__":
    main_cli()
