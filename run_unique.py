import sys
import glob
from AaronTools.job_control import SubmitProcess

def run_jobs(spec = False):
    if spec == True:
        com_files = glob.glob("*.com") + glob.glob("*.in") + glob.glob("*.inp")
    else:
        com_files = glob.glob("unique_files/*.com") + glob.glob("unique_files/*.in") + glob.glob("unique_files/*.inp")
    if not com_files:
        print("No runnable files found in current directory.")
        sys.exit(1)

    for f in com_files:
        submit_process = SubmitProcess(f, 12, 8, 12)

         # submit job
        try:
            submit_process.submit()
        except Exception as e:
            print(f"failed to submit {f}: {e}")
 
