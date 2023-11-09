import subprocess
import os
import shutil
import time


def timed_function(func, *args, **kwargs):
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    print(f"Time taken for {func.__name__}: {end_time - start_time:.2f} seconds")
    return result


# Get the current directory of the test script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define the relative path to your Python script
MYSCRIPT = "PlastomeOldCode.py"

# Combine the directory of the test script with the relative path to the Python script
full_script_path = os.path.join(script_dir, MYSCRIPT)

# Define the path to the output directory
folder_INT = "benchmarking1/output_INT"

# Step 1: Extract the tar.gz file
timed_function(subprocess.run, ["tar", "-xvf", "benchmarking1.tar.gz"])

# Step 2: Change the current directory to benchmarking1/
os.chdir("benchmarking1")

# Step 3: Create necessary directories
subprocess.run(["mkdir", "-p", folder_INT])
subprocess.run(["mkdir", "-p", f"{folder_INT}/01_unalign"])
subprocess.run(["mkdir", "-p", f"{folder_INT}/02_aligned"])
subprocess.run(["mkdir", "-p", f"{folder_INT}/02_aligned/fasta"])
subprocess.run(["mkdir", "-p", f"{folder_INT}/02_aligned/nexus"])

# Step 4: Run your Python script using the full_script_path
timed_function(
    subprocess.run,
    ["python", full_script_path, "-i", ".", "-o", folder_INT, "-s", "int"],
)

# run this to remove the folder, if not can comment out
# Step 5: Delete the benchmarking1 directory and its contents
os.chdir(script_dir)  # Move back to the directory containing the test script
shutil.rmtree("benchmarking1")
print("benchmarking1 directory and its contents have been deleted.")
