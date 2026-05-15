import json
import os
import sys
import shutil
import subprocess
from datetime import datetime

def load_config(config_path="config.json"):
    """reads the JSON file and returns it as a dictionary"""
    if not os.path.exists(config_path):
        print(f"Error: Could not find {config_path}")
        sys.exit(1)

    with open(config_path, "r") as file:
        config = json.load(file)

    return config

def create_run_environment(cfg):
    """creates a unique directory for the run and copies the data"""

    timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
    run_index = 1
    run_dir = f"runs/{cfg['target_name'].lower()}_{run_index}"
    while os.path.isdir(run_dir):
        run_index += 1
        run_dir = f"runs/{cfg['target_name'].lower()}_{run_index}"
    # latest_symlink = "runs/latest"
    # if os.path.islink(latest_symlink):
    #     os.unlink(latest_symlink)

    # os.symlink(os.path.basename(run_dir), latest_symlink)
    template_dir = "templates"

    copy_files = ["Data.h",
"Data.cpp",
f"{cfg['model_type']}.h",
f"{cfg['model_type']}.cpp",
"main.cpp",
"Makefile"
    ]

    # make directory
    os.makedirs(run_dir, exist_ok=True)
    print(f"Making run directory targeting {cfg['target_name']} using data {cfg['data_file']}")
    print(f"Created run directory: {run_dir}")

    # copy data file into run dir
    print("Starting file copies")
    data_source = cfg["data_file"]
    if os.path.exists(data_source):
        shutil.copy(data_source, os.path.join(run_dir, "data.txt"))
        print(f"Copied data from {data_source} to {run_dir}/data.txt")
    else:
        print(f"Data file {data_source} not found")

    if os.path.exists("config.json"):
        shutil.copy("config.json", os.path.join(run_dir, f"{timestamp}_config.json"))
        print(f"Successfully copied config file")
    else:
        print(f"Config file not found")

    for filename in copy_files:
        source_path = os.path.join(template_dir, filename)
        dest_path = os.path.join(run_dir, filename)

        if os.path.exists(source_path):
            shutil.copy(source_path, dest_path)
            print(f"Copied file {filename} from {source_path} to {dest_path}")
        else:
            print(f"File {filename} not found")

    print("Compiling model")
    try:
        subprocess.run(["make"], cwd=run_dir, check=True)
        print("Compilation success")
    except subprocess.CalledProcessError:
        print("Compilation failed")

    return run_dir

def generate_options_file(run_dir, cfg):
    """automatically writes json config into DNest4 OPTIONS file"""
    options_path = os.path.join(run_dir, "OPTIONS")

    options_content = f"""
    {cfg["num_particles"]}
    {cfg["new_level_interval"]}
    {cfg["save_interval"]}
    {cfg["thread_steps"]}
    {cfg["max_num_levels"]}
    {cfg["lambda"]}
    {cfg["beta"]}
    {cfg["max_num_saves"]}
    """

    with open(options_path, "w") as f:
        f.write(options_content)

    print("Generated OPTIONS file")

if __name__ == "__main__":
    cfg = load_config()
    print("=============================")

    # run generator
    current_run_dir = create_run_environment(cfg)
    generate_options_file(current_run_dir, cfg)
    print(f"Success! New run created at {current_run_dir}")

    print("=============================")