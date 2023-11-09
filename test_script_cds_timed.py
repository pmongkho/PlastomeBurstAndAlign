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
MYSCRIPT = "PlastomeRegionBurstAndAlign.py"

# Combine the directory of the test script with the relative path to the Python script
full_script_path = os.path.join(script_dir, MYSCRIPT)

# Define the path to the output directory
folder_CDS = "benchmarking1/output_CDS"

# Wrap each step that you want to time with the timed_function

# Step 1: Extract the tar.gz file
timed_function(subprocess.run, ["tar", "-xvf", "benchmarking1.tar.gz"])

# Step 2: Change the current directory to benchmarking1/
os.chdir("benchmarking1")

# Step 3: Create necessary directories
timed_function(subprocess.run, ["mkdir", "-p", folder_CDS])
timed_function(subprocess.run, ["mkdir", "-p", f"{folder_CDS}/01_unalign"])
timed_function(subprocess.run, ["mkdir", "-p", f"{folder_CDS}/02_aligned"])
timed_function(subprocess.run, ["mkdir", "-p", f"{folder_CDS}/02_aligned/fasta"])
timed_function(subprocess.run, ["mkdir", "-p", f"{folder_CDS}/02_aligned/nexus"])

# Step 4: Run your Python script using the full_script_path
timed_function(
    subprocess.run,
    ["python", full_script_path, "-i", ".", "-o", folder_CDS, "-s", "cds"],
)

# Step 5: Delete the benchmarking1 directory and its contents
# Make sure you really want to delete before uncommenting the next lines
os.chdir(script_dir)  # Move back to the directory containing the test script
timed_function(shutil.rmtree, "benchmarking1")
print("benchmarking1 directory and its contents have been deleted.")
