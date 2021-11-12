import cclib
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# gaussian_scan_output = "/home/unibas/boittier/RDKit_G2/amide1.pdb/SCAN_1_2_3_4_S_36_10.0/SCAN_amide1.pdb.out"
# data_path = "/data/unibas/boittier/fdcm_test/"
# csv_out_name = "fdcm_test.csv"

gaussian_scan_output = "/home/unibas/boittier/RDKit_G2/amide1.pdb/SCAN_1_2_3_4_S_36_10.0/SCAN_amide2.pdb.out"
data_path = "/data/unibas/boittier/amide2/"
csv_out_name = "amide2.csv"

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

frames = []
errors = []
a1_ = []
a2_ = []
d1_ = []
for f in files:
    print(f)
    result = open(f).readlines()[-1].split()[1]
    frame = int(f.split("/")[-2].split("_")[1])
    error = float(result)
    frames.append(frame)
    errors.append(error)
    frame = frame -1
    a1_.append(a1[frame])
    a2_.append(a2[frame])
    d1_.append(d1[frame])


print(frames, len(a1_), len(a2_), len(d1_), len(frames))
df = pd.DataFrame({"frame": frames, "error": errors, "a1":a1_, "a2": a2_, "d1": d1_})


print(df)
df.to_csv(csv_out_name)



# df.plot("frame", "error", kind="bar")
# plt.savefig("test.png")


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# a1_ = np.array(a1_).reshape(len(a1_))
# a2_ = np.array(a2_).reshape(len(a1_))
# d1_ = np.array(d1_).reshape(len(a1_))

# ax.scatter(a1_, a2_, d1_, c=errors)

# sa1 = list(set(a1_))
# sa1.sort()
# sa2 = list(set(a2_))
# sa2.sort()
# sd1 = list(set(d1_))
# sd1.sort()


# texts = []
# for i, _ in enumerate(a1_):
#     text = f"{sa1.index(a1_[i])}{sa2.index(a2_[i])}{sd1.index(d1_[i])}"
#     texts.append(text)
#     ax.text(a1_[i], a2_[i], d1_[i], text, "y")

# ax.set_xlabel("a1")
# ax.set_ylabel("a2")
# ax.set_zlabel("d1")


# plt.savefig("3d.png")
# plt.clf()

# plt.bar(texts, errors)
# plt.savefig("bar2.png")


