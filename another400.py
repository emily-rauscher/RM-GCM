#!/usr/bin/env python

import os
import pathlib
import shutil
import glob
import numpy as np
NDAYS = 200
# FOLDERS = ["Planet_Run","Planet_Run_Restart","Planet_Run_Restart2","Planet_Run_Restart3","Planet_Run_Restart4","Planet_Run_Restart5","Planet_Run_Restart6","Planet_Run_Restart7","Planet_Run_Restart8","Planet_Run_Restart9","Planet_Run_Restart10"]
# i=len(FOLDERS)-1
# while not os.path.exists(FOLDERS[i]+"/fort.29"):
#     i-=1

FOLDER1 = os.getcwd()
# FOLDER2 = FOLDERS[i+1]

#shutil.copy(FOLDER1+"/fort.7", "./")

with open(FOLDER1+"/fort.29") as f:
    lines = f.readlines()
    KOUNT = float(lines[-1].split()[-1])
with open('fort.7') as f:
    lines=f.readlines()
for line in lines:
    if "KOUNTR" in line:
        KOUNTR = float(line.split()[-1][:-1])
#         print(line)
    elif "TSPD" in line:
        TSPD = float(line.split()[-1][:-1])
#         print(line)
    elif "KRUN" in line:
        KRUN = float(line.split()[-1][:-1])
#         print(line)
    elif "BEGDAY" in line:
        BEGDAY = float(line.split()[-1][:-1])
KRUNNEW = NDAYS * TSPD
BEGDAYNEW = np.round(KOUNT / TSPD)

#print(BEGDAY)
print(KRUNNEW, BEGDAYNEW, TSPD, FOLDER1)


