import pandas as pd
from ase import Atoms
from dscribe.descriptors import SOAP
import numpy as np
import sys
from ARS import *


def do_transformation(old_global_xyz, cube, frame_file):
    ars_obj = ARS(old_global_xyz, cube, frame_file, pcube_2=cube)

    f = open("frame_33_local.xyz").readlines()
    fitted_local = []
    for x in f[2:]:
        fitted_local.append([float(x.split()[1]), float(x.split()[2]), float(x.split()[3])])

    ars_obj.set_local_charge_positions(np.array(fitted_local))
    cp = ars_obj.local_to_global()
    ars_obj.set_charge_positions_plus(cp)
    ars_obj.save_charges_global("kernel.xyz")


BOHR_TO_ANGSTROM = 0.529177
model = pd.read_pickle(r"best_model.pkl")


def z_to_char(z):
    if z == 1:
        return "H"
    elif z == 6:
        return "C"
    elif z == 8:
        return "O"


header = """24
s                      x[A]                      y[A]                      z[A]                      q[e]\n"""
line_format = "N {} {} {} {} \n"

#  Hard coding charges for ester
q = [0.6959175615854226, 0.4294665048044459, 0.7289443612696791, -0.8899528116456963, .4821174903544280,
     0.6035297564939646, 0.0604057803354635, 6991174418445399, 0.4131324281244263, .9286014603942530,
     -0.1057946580878595, 9961038253671308, .3737552305035935, .9714181100989372, 0.9973398242463560,
     0.7637225416907655, -0.9974375214017446, 0.4871145088925188, .4662300369825688, .6223673236702755,
     0.9509229066582735, -0.9922381659518612, -0.6086911331096487, 0.3897866237988844]


def cube_to_xyz(path):
    lines = open(path).readlines()
    n_atoms = int(lines[2].split()[0])
    atoms = []
    for line in lines[7 - 1:7 + n_atoms - 1]:
        Z, x, y, z = line.split()[1:]
        atoms.append([float(Z), float(x) * BOHR_TO_ANGSTROM, float(y) * BOHR_TO_ANGSTROM, float(z) * BOHR_TO_ANGSTROM])
    return atoms


path = "/home/boittier/Documents/PhD/data/ester/frame_33.chk.d.cube"


def model_to_local(cube_path, local_output):
    atoms = cube_to_xyz(cube_path)
    print(atoms)
    RCUT = 5
    NMAX = 4
    LMAX = 5
    #  SOAP descriptors
    soap_desc = SOAP(species=["C", "H", "O"], rcut=RCUT, nmax=NMAX, lmax=LMAX, crossover=True)
    atom_types = [z_to_char(x[0]) for x in atoms]
    string = [x[1:] for x in atoms]
    atoms = Atoms(atom_types, string)
    samples = [atoms]
    der, des = soap_desc.derivatives(samples, method="analytical", return_descriptor=True, n_jobs=8)
    print(des.shape)
    a_, b_, c_ = des.shape
    X = des[0].reshape((1, b_ * c_))
    pred = model.predict(X)
    pred = pred[0].reshape(24, 3)
    f = open(local_output, "w")
    f.write(header)
    for i, xyz in enumerate(pred):
        f.write(line_format.format(*xyz, q[i]))


        
frame_file = "/home/boittier/ester_traj/frames.txt"
cube = "/home/boittier/ester_traj/t0/frame_33.chk.p.cube"
old_global_xyz = "/home/boittier/FDCM/ester_t1_100/frame_33/refined.xyz"
local_output = "frame_33_local.xyz"

def main():
    frame_file = sys.argv[1]
    cube = sys.argv[2]
    old_global_xyz = sys.argv[3]
    local_output = sys.argv[4]
    model_to_local(path, local_output)
    do_transformation(old_global_xyz, cube, frame_file)

if __name__ == "__main__":
    main()