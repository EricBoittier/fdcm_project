import os

import cclib
import numpy as np
import pandas as pd


def get_neighbours(point, ranges):
    neighbours = []
    for i, dim in enumerate(point):
        if ranges[i][0] <= (dim - 1) <= ranges[i][1]:
            n = point.copy()
            n[i] = dim - 1
            neighbours.append(tuple(n))
        if ranges[i][0] <= (dim + 1) <= ranges[i][1]:
            n = point.copy()
            n[i] = dim + 1
            neighbours.append(tuple(n))

    return neighbours


def add_key_int(morton):
    s_a1 = list(set(morton["a1"]))
    s_a1.sort()
    s_a2 = list(set(morton["a2"]))
    s_a2.sort()
    s_d1 = list(set(morton["d1"]))
    s_d1.sort()
    morton["a1_"] = morton["a1"].apply(s_a1.index)
    morton["a2_"] = morton["a2"].apply(s_a2.index)
    morton["d1_"] = morton["d1"].apply(s_d1.index)
    return morton


def key_to_frame(df, i):
    _ = df[df["a1_"] == i[0]]
    _ = _[_["a2_"] == i[1]]
    _ = _[_["d1_"] == i[2]]
    return int(list(_["frame"])[0])


def get_local_charges(path):
    tmp_dict = {}
    local_file = open(path).readlines()
    n_lines = int(local_file[0])
    charges = local_file[2:n_lines + 2]
    charges = np.array([[float(x.split()[1]),
                         float(x.split()[2]),
                         float(x.split()[3]),
                         float(x.split()[4])] for x in charges])

    for i, axis in enumerate(["x", "y", "z", "q"]):
        for c in range(n_lines):
            if f"{axis}_c{c}" not in tmp_dict.keys():
                tmp_dict[f"{axis}_c{c}"] = charges[c, i]

    return tmp_dict


def get_path_neighbours(args):
    gaussian_scan_output = args.gaussian_scan_output
    print(gaussian_scan_output)
    p = cclib.io.ccopen(gaussian_scan_output)
    p = p.parse()
    scan_names = p.scannames
    scan_parms = np.array(p.scanparm)
    print(scan_parms)

    a1 = [x // 1 for x in scan_parms[0]]
    a2 = [x // 1 for x in scan_parms[1]]
    d1 = [x // 1 for x in scan_parms[2]]

    a1 = {i: x for i, x in enumerate(a1)}
    a2 = {i: x for i, x in enumerate(a2)}
    d1 = {i: x for i, x in enumerate(d1)}

    a1_ = []
    a2_ = []
    d1_ = []
    frames = []

    for frame in range(len(a1)):
        frames.append(frame)
        a1_.append(a1[frame])
        a2_.append(a2[frame])
        d1_.append(d1[frame])

    df = pd.DataFrame({"frame": frames, "a1": a1_, "a2": a2_, "d1": d1_})
    df = add_key_int(df)
    ranges = [[0, len(set(a1_))], [0, len(set(a2))], [0, len(set(d1))]]

    path = []
    neighbours = []

    for key, j in enumerate(range(len(frames))):
        r = df[df["frame"] == key]
        i = np.array([int(r["a1_"]), int(r["a2_"]), int(r["d1_"])])

        ns = get_neighbours(i, ranges)
        frame = key_to_frame(df, i)
        print("visiting: frame_", frame)
        print(ns)
        # visited.append(frame)

        n = [key_to_frame(df, ii) for ii in ns]
        neighbours.append(n)
        path.append(key)


def analyse(args):
    data_path = args.output_dir
    csv_out_name = args.csv_out_name

    a1 = [x // 1 for x in scan_parms[0]]
    a2 = [x // 1 for x in scan_parms[1]]
    d1 = [x // 1 for x in scan_parms[2]]

    a1 = {i: x for i, x in enumerate(a1)}
    a2 = {i: x for i, x in enumerate(a2)}
    d1 = {i: x for i, x in enumerate(d1)}

    files = [os.path.join(data_path, x, "GD.log") for x in os.listdir(data_path) if x.__contains__("frame")]
    local_file_dirs = [os.path.join(data_path, x) for x in os.listdir(data_path) if
                   x.__contains__("frame")]
    local_files = []
    for lfdir in local_file_dirs:
        local_file = [x for x in os.listdir(lfdir) if x.endswith(".local")][0]
        local_files.append(os.path.join(lfdir, local_file))

    frames = []
    errors = []
    local_charges = []

    for f, local in zip(files, local_files):
        local_charges.append(get_local_charges(local))
        result = open(f).readlines()[-1].split()[1]
        frame = int(f.split("/")[-2].split("_")[1])
        error = float(result)
        frames.append(frame)
        errors.append(error)

    lc_df = pd.DataFrame(local_charges)

    df = pd.DataFrame({"frame": frames, "error": errors})
    df = df.join(lc_df)
    print(df)
    df.to_csv(csv_out_name)
    # df = add_key_int(df)



def analyse_gaussian(args):
    gaussian_scan_output = args.gaussian_scan_output
    data_path = args.data_path
    csv_out_name = args.csv_out_name
    print(gaussian_scan_output)
    p = cclib.io.ccopen(gaussian_scan_output)
    p = p.parse()
    scan_names = p.scannames
    scan_parms = np.array(p.scanparm)
    print(scan_parms)

    a1 = [x // 1 for x in scan_parms[0]]
    a2 = [x // 1 for x in scan_parms[1]]
    d1 = [x // 1 for x in scan_parms[2]]

    a1 = {i: x for i, x in enumerate(a1)}
    a2 = {i: x for i, x in enumerate(a2)}
    d1 = {i: x for i, x in enumerate(d1)}

    files = [os.path.join(data_path, x, "GD.log") for x in os.listdir(data_path) if x.__contains__("frame")]
    local_files = [os.path.join(data_path, x, "refined.xyz.local") for x in os.listdir(data_path) if
                   x.__contains__("frame")]

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
    df = pd.DataFrame({"frame": frames, "error": errors, "a1": a1_, "a2": a2_, "d1": d1_})
    df = df.join(lc_df)
    print(df)
    df.to_csv(csv_out_name)

    df = add_key_int(df)