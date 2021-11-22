import pandas as pd
from ase import Atoms
from dscribe.descriptors import SOAP
import numpy as np
import sys
from ARS import *


def func(data, a, b,c,d,e,f,g,h,i,j):
    return a*np.sin(b*data[0]+c)+ d + e*data[1]**3 + f*data[1]**2 + g*data[1] + h*data[2]**3 + i*data[2]**2 + j*data[2]

def func2(data, a, b,c,d,g,j):
    return a*np.sin(b*data[0]+c)+ d + g*data[1]**2 + j*data[2]**2

def func3(data, a, b,c,d):
    return a*np.sin(b*data[0]+c)+ d #+ g*data[1] + j*data[2]


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
model = pd.read_pickle(r"kernel_data/best_model.pkl")


def z_to_char(z):
    if z == 1:
        return "H"
    elif z == 6:
        return "C"
    elif z == 7:
        return "N"
    elif z == 8:
        return "O"


header = """24
s                      x[A]                      y[A]                      z[A]                      q[e]\n"""
line_format = "N {} {} {} {} \n"


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

        
def fit_to_local(cube_path, fit_path, old_mdcm_path, local_output, frame_file):
    #  Load the fit
    per_charge_params = pd.read_pickle(fit_path)
    
    #  Creat an ARS object for atom positions, charges, angles, dih, etc
    ars_obj = ARS(old_mdcm_path, cube_path, frame_file, pcube_2=cube_path)
    
    q = ars_obj.c_charges
    n_charges = len(q)
    a1 = ars_obj.get_angle(0, 1, 7)
    print(a1)
    a2 = ars_obj.get_angle(1, 7, 9)
    print(a2)
    d1 = ars_obj.get_dih( 0, 1, 7, 9)
    print(d1)
    
    #  Calculate the parameters from the cube file
    f = open(local_output, "w")
    f.write(header)
    fitted_local = []
    for i in range(n_charges):
        x = func([d1, a1, a2], *per_charge_params[i])
        y = func([d1, a1, a2], *per_charge_params[i+n_charges])
        z = func([d1, a1, a2], *per_charge_params[i+n_charges*2])
        f.write(line_format.format(x, y, z, q[i]))
        fitted_local.append([x, y, z])
        
    ars_obj.set_local_charge_positions(np.array(fitted_local))
    cp = ars_obj.local_to_global()
    ars_obj.set_charge_positions_plus(cp)
    ars_obj.save_charges_global("fit.xyz")

    pass

        
# frame_file = "/home/boittier/ester_traj/frames.txt"
# cube = "/home/boittier/ester_traj/t0/frame_33.chk.p.cube"
# old_global_xyz = "/home/boittier/FDCM/ester_t1_100/frame_33/refined.xyz"
# local_output = "frame_33_local.xyz"

def main():
#     frame_file = sys.argv[1]
#     cube = sys.argv[2]
#     old_global_xyz = sys.argv[3]
#     local_output = sys.argv[4]

#     model_to_local(cube, local_output)
#     do_transformation(old_global_xyz, cube, frame_file)
    frame_file = "/home/unibas/boittier/fdcm_project/mdcms/amide/model1/frames.txt"
    cube_path = "/data/unibas/boittier/fdcm/amide/scan-large/SCAN_amide1.pdb-0.xyz.chk.d.cube"
    local_output = "test_fit.xyz"
    fit_path = "/home/unibas/boittier/fdcm_project/fdcm_notebooks/loop-model.pkl"
    old_mdcm_path = "/home/unibas/boittier/fdcm_project/mdcms/amide/model1/24_charges_refined.xyz"
    
    fit_to_local(cube_path, fit_path, old_mdcm_path, local_output, frame_file)
    

if __name__ == "__main__":
    main()