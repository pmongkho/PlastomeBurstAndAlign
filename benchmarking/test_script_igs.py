import subprocess
import os
import shutil
from security import safe_command

# Get the current directory of the test script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define the relative path to your Python script
MYSCRIPT = "../PlastomeRegionBurstAndAlign.py"

# Combine the directory of the test script with the relative path to the Python script
full_script_path = os.path.join(script_dir, MYSCRIPT)

# Define the path to the output directory
folder_IGS = "benchmarking1/output_IGS"

# Step 1: Extract the tar.gz file
safe_command.run(subprocess.run, ["tar", "-xvf", "benchmarking1.tar.gz"])

# Step 2: Change the current directory to benchmarking1/
os.chdir("benchmarking1")

# Step 3: Create necessary directories
safe_command.run(subprocess.run, ["mkdir", "-p", folder_IGS])
safe_command.run(subprocess.run, ["mkdir", "-p", f"{folder_IGS}/01_unalign"])
safe_command.run(subprocess.run, ["mkdir", "-p", f"{folder_IGS}/02_aligned"])
safe_command.run(subprocess.run, ["mkdir", "-p", f"{folder_IGS}/02_aligned/fasta"])
safe_command.run(subprocess.run, ["mkdir", "-p", f"{folder_IGS}/02_aligned/nexus"])

# Step 4: Run your Python script using the full_script_path
safe_command.run(subprocess.run, ["python", full_script_path, "-i", ".", "-o", folder_IGS, "-s", "igs", "-t", "5", "-l", "6"])

# run this to remove the folder, if not can comment out
# Step 5: Delete the benchmarking1 directory and its contents
os.chdir(script_dir)  # Move back to the directory containing the test script
shutil.rmtree("benchmarking1")
print("benchmarking1 directory and its contents have been deleted.")

