#!/usr/bin/env python3

import os
from hfrpkg.multi.fill import fill_logs
from hfrpkg.multi.comp import main as comp

def combined_postprocess():
    
    fill_logs("unique_files")

    comp()

    print("Done.")

def main_cli():
    combined_postprocess()


if __name__ == "__main__":
    main_cli()
