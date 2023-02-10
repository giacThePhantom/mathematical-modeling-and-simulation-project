import os
import sys
import glob

file_list = glob.glob("../params_files/*.csv")
for file in file_list:
    os.remove(file)

min = int(sys.argv[1])
step = int(sys.argv[2])
max = int(sys.argv[3])
iterations = int((max-min)/step) + 2

k = min
for i in range(0, iterations):
    ca = min
    for y in range(0, iterations):
        new_file = open(f"../params_files/{int(k)}_{int(ca)}.csv", "w")
        cols = "tmax,nktot,mcatot,phi_m,phi_n,va,vb,vc,vd,vk,vl,vca,gk,gl,gca,c,v0,m0,n0,i\n"
        params = [2000.0,float(k),float(ca),0.4,0.04,-1.2,18.0,2.0,30.0,-84.0,-60.0,120.0,8.0,2.0,4.4,20.0,-50.0,0.0,int(float(k)/2),100.0]
        new_file.writelines([cols, str(params)[1:len(str(params))-1].replace(" ", "")])
        new_file.close()
        if (y == 0):
            if (min > step):
                ca = min
            else:
                ca = step
        else:
            ca += step
    if (i == 0):
        if (min > step):
            k = min
        else:
            k = step
    else:
        k += step