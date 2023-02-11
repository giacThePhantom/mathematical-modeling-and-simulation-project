import re
import os
import sys
import glob
import subprocess

from tqdm import tqdm
from pathlib import Path
from itertools import combinations

# Input arguments
# argv[1] rtc or gill
# argv[2] name of the log file
# argv[3] max parallel jobs

test_list = [os.path.join("params_files", x) for x in os.listdir("../params_files")]
pairs_list = list(combinations(test_list, r=2))[0:10] 
processes = []
j = 0

function = {
    "rtc": "rtc_crp",
    "gill": "gill_crn"
}

f = open(os.path.join("..", "results", "distance", f"log_{sys.argv[2]}.txt"), "a")
f.write("sim_name,result\n")

s_p = 0;
for files in pairs_list:

    # Create log file for this run
    f_file = re.search(r"[0-9]{1,}_[0-9]{1,}", files[0]).group(0)
    s_file = re.search(r"[0-9]{1,}_[0-9]{1,}", files[1]).group(0)
    sim_name = f"{f_file}-with-{s_file}"

    # Pass correct file paths to matlab
    file1_path = os.path.join("..", "comparisons", files[0])
    file2_path = os.path.join("..", "comparisons", files[1])

    # Run matlab process 
    cmd = f"""matlab -nojvm -batch "{function[sys.argv[1]]}('{file1_path}','{file2_path}')" """
    p = subprocess.Popen(cmd, 
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        cwd="../../src/")
    processes.append((p, sim_name))
    
    # Set how many jobs are allowed to run at the same time
    j += 1
    if j % int(sys.argv[3]) == 0:
        print(f"Processing {s_p} to {j}...")
        p.wait()
        s_p = j


for p, sim_name in processes:
    p.wait()
    out, err = p.communicate()
    if (not(err)):
        decoded_out = out.decode("utf-8")
        f.write(f"{sim_name},{decoded_out}\n")
    else:
        decoded_ett = err.decode("utf-8")
        f.write(f"{sim_name},{decoded_ett}\n")
f.close()
