import sys
import subprocesses
import os

processes = []
j = 0
for i in arg_list:
    cmd = f'matlab -r "{sys.argv[1]} {i[1]} {i[2]}"'
    f = os.tmpfile()
    p = subprocess.Popen(cmd, stdout = f)
    processes.append((p,f))
    j += 1
    if j % 10 == 0:
        p.wait()

for p, f in processes:
    p.wait()
    f.seek(0)
    print(f.read())
