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

def get_process_output(p, sim_name, f):
        out, err = p.communicate()
        if (not(err)):
            decoded_out = out.decode("utf-8")
            print(f"{sim_name},{decoded_out}\n")
            f.write(f"{sim_name},{decoded_out}\n")
        else:
            decoded_err = err.decode("utf-8")
            print(f"{sim_name},{decoded_err}\n")
            f.write(f"{sim_name},{decoded_err}\n")


test_list = [os.path.join("params_files", x) for x in os.listdir("../params_files")]
pairs_list = list(combinations(test_list, r=2))[0:12] 
processes = []
j = 0

function = {
    "rtc": "rtc_crp",
    "gill": "gill_crn"
}

f = open(os.path.join("..", "results", "distance", f"log_{sys.argv[2]}.txt"), "a")
f.write("sim_name,result\n")

s_p = 0;
done = []
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
        print(f"Processing {s_p} to {j}..")
        for i in range(s_p, j):
            get_process_output(p, sim_name, f)
            done.append(i)
            s_p = j

#Left over
for i in range(0, len(processes)):
    if (i not in done):
        p, sim_name = processes[i]
        print(f"Processing {sim_name}..")
        get_process_output(p, sim_name, f)
f.close()
