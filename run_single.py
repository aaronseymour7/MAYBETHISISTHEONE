import sys
import glob
from AaronTools.job_control import SubmitProcess

def run_jobs():
    def get_extension(index_path="index.txt"):
        ext_map = {
            "gaussian": ".com",
            "orca": ".inp",
            "psi4": ".in"
        }
        try:
            with open(index_path) as f:
                first_line = f.readline()
                # Expecting format: "Level:\t reaction_fn_name\t software: \t software_name"
                parts = first_line.strip().split("\t")
                if len(parts) >= 4:
                    return ext_map.get(parts[3].lower(), None)
        except FileNotFoundError:
            pass
        return None
    ext = get_extension()
    if ext is None:
        print("Could not determine file extension from index.txt")
        sys.exit(1)
    com_files = glob.glob("*"+ ext)
    if not com_files:
        print(f"No {ext} files found in current directory.")
        sys.exit(1)

    for f in com_files:
        submit_process = SubmitProcess(f, 12, 8, 12)

         # submit job
        try:
            submit_process.submit()
        except Exception as e:
            print(f"failed to submit {f}: {e}")
def main_cli():
    run_jobs()

if __name__ == "__main__":
    main_cli()
