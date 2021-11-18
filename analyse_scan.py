import cclib
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_local_charges(path):
    tmp_dict = {}
    local_file = open(path).readlines()
    n_lines = int(local_file[0])
    charges = local_file[2:n_lines+2]
    charges = np.array([ [float(x.split()[1]),
              float(x.split()[2]),
              float(x.split()[3]),
              float(x.split()[4])] for x in charges])

    for i, axis in enumerate(["x", "y", "z", "q"]):
        for c in range(n_lines):
            if f"{axis}_c{c}" not in tmp_dict.keys():
                tmp_dict[f"{axis}_c{c}"] = charges[c,i]
            
            
    return tmp_dict


# gaussian_scan_output = "/home/unibas/boittier/RDKit_G2/amide1.pdb/SCAN_1_2_3_4_S_36_10.0/SCAN_amide1.pdb.out"
# data_path = "/data/unibas/boittier/morton_test"
# csv_out_name = "morton_test.csv"

gaussian_scan_output = "/home/unibas/boittier/RDKit_G2/amide1.pdb/SCAN_1_2_3_4_S_36_10.0/test/SCAN_amide1.pdb.out"
data_path = "/data/unibas/boittier/morton-large-control5"
csv_out_name = "morton-large-control5.csv"

#gaussian_scan_output = "/home/unibas/boittier/RDKit_G2/amide1.pdb/SCAN_1_2_3_4_S_36_10.0/SCAN_amide2.pdb.out"
#data_path = "/data/unibas/boittier/amide2/"
#csv_out_name = "amide2.csv"

p = cclib.io.ccopen(gaussian_scan_output)
p = p.parse()
scan_names = p.scannames
scan_parms = np.array(p.scanparm)
print(scan_parms)

a1 = [x//1 for x in scan_parms[0]]
a2 = [x//1 for x in scan_parms[1]]
d1 = [x//1 for x in scan_parms[2]]

a1 = {i : x for i, x in enumerate(a1)}
a2 = {i : x for i, x in enumerate(a2)}
d1 = {i : x for i, x in enumerate(d1)}


files = [os.path.join(data_path, x, "GD.log") for x in os.listdir(data_path) if x.__contains__("frame")]
local_files =  [os.path.join(data_path, x, "refined.xyz.local") for x in os.listdir(data_path) if x.__contains__("frame")]

frames = []
errors = []
a1_ = []
a2_ = []
d1_ = []
local_charges = []

for f, local in zip(files, local_files):
    local_charges.append(get_local_charges(local))
    result = open(f).readlines()[-1].split()[1]
    frame = int(f.split("/")[-2].split("_")[1])
    error = float(result)
    frames.append(frame)
    errors.append(error)
    a1_.append(a1[frame])
    a2_.append(a2[frame])
    d1_.append(d1[frame])

lc_df = pd.DataFrame(local_charges)
print(lc_df)
print(frames, len(a1_), len(a2_), len(d1_), len(frames))
df = pd.DataFrame({"frame": frames, "error": errors, "a1":a1_, "a2": a2_, "d1": d1_})
df = df.join(lc_df)
print(df)
df.to_csv(csv_out_name)



