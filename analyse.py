import sys
import os
import dnest4.classic as dn4

def analyse(run_dir):
    original_dir = os.getcwd()

    try:
        os.chdir(f"runs/{run_dir}")
        print(f"=== Analysing {run_dir} ===")

        dn4.postprocess()

        print(f"=== Analysis complete ===")
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 analyse.py <path_to_run_directory>")
        sys.exit(1)

    run_directory = sys.argv[1]
    analyse(run_directory)