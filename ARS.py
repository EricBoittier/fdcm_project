import argparse

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from scipy.spatial import distance
from scipy.spatial.transform import Rotation as Kabsch

from Cube import read_cube

BOHR_TO_ANGSTROM = 0.529177


def angle(a, b, c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    a3 = np.degrees(np.arccos(cosine_angle))
    return a3


# This is the straightforward approach as outlined in the answers to
# "How do I calculate a dihedral angle given Cartesian coordinates?"
def dihedral2(p):
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array([v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1, 1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2(y, x))


def read_cube_file(filepath):
    pcube_data, pcube_meta = read_cube(filepath)
    ap = []
    an = []
    for i in pcube_meta["atoms"]:
        atom = list(i[1])
        ap.append([x * BOHR_TO_ANGSTROM for x in atom[1:]])
        an.append(atom[1])
    return ap, an


def read_mdcm_xyz(filepath):
    xyz_file = open(filepath).readlines()
    n_charges = int(xyz_file[0])
    #  read number of charges from first line (xyz format)
    charge_lines = xyz_file[2:n_charges + 2]
    # Read atoms and charges
    c_positions = []
    c_charges = []
    for charge in charge_lines:
        on, x, y, z, c = charge.split()
        c_positions.append([float(x), float(y), float(z)])
        c_charges.append(float(c))
    return c_positions, c_charges


def get_local_axis(atom_pos, frame_atoms, method="bond"):
    """
    method: "bond" z-axis a-b
            "bisector" z-axis = bisector of a-b,b-c
    Inputs:
                atom_positions, frames
    Returns:
                List of Lists of Frame Vectors [ [x_v, y_v, z_v], ...  ] in order of frames
    """
    n_frames = len(frame_atoms)
    frame_vectors = []
    for f in range(n_frames):
        a_index, b_index, c_index = frame_atoms[f]
        a, b, c = frame_atoms[f]
        # adjust indexing
        a = atom_pos[a - 1]
        b = atom_pos[b - 1]
        c = atom_pos[c - 1]
        distance_ab = distance.euclidean(a, b)
        b1_x = (a[0] - b[0]) / distance_ab
        b1_y = (a[1] - b[1]) / distance_ab
        b1_z = (a[2] - b[2]) / distance_ab

        distance_bc = distance.euclidean(b, c)
        b2_x = (b[0] - c[0]) / distance_bc
        b2_y = (b[1] - c[1]) / distance_bc
        b2_z = (b[2] - c[2]) / distance_bc

        #  Z axes
        ez1 = np.array([b1_x, b1_y, b1_z])
        ez3 = np.array([b2_x, b2_y, b2_z])

        if method == "bond":
            ez2 = np.array([b1_x, b1_y, b1_z])

        elif method == "bisector":
            """ Calculate Z(2) as bisector
            """
            bi_x = ez1[0] + ez3[0]
            bi_y = ez1[1] + ez3[1]
            bi_z = ez1[2] + ez3[2]

            #  get norm
            r_bi = np.sqrt(bi_x ** 2 + bi_y ** 2 + bi_z ** 2)
            #  normalize
            ez2 = np.array([bi_x, bi_y, bi_z]) / r_bi

            if r_bi < 0.0001:
                print("Colinearity detected! (Bad)")

        else:
            assert False, "No valid method supplied!"

        #  Y axes
        ey1 = np.zeros(3)
        ey1[0] = b1_y * b2_z - b1_z * b2_y
        ey1[1] = b1_z * b2_x - b1_x * b2_z
        ey1[2] = b1_x * b2_y - b1_y * b2_x
        re_x = np.sqrt(ey1[0] ** 2 + ey1[1] ** 2 + ey1[2] ** 2)
        ey1[0] = ey1[0] / re_x
        ey1[1] = ey1[1] / re_x
        ey1[2] = ey1[2] / re_x

        #         ey1 = np.cross(ez1, ez3)

        ey2 = ey1
        ey3 = ey1

        #  X axes
        ex1 = np.zeros(3)
        ex3 = np.zeros(3)
        #  ex1 and ex2
        ex1 = np.cross(ey1, ez1)
        if method == "bond":
            ex2 = ex1
        else:
            ex2 = np.cross(ey2, ez2)
        #  ex3
        ex3 = np.cross(ey3, ez3)

        frame_vectors.append(([ex1, ey1, ez1],
                              [ex2, ey2, ez2],
                              [ex3, ey3, ez3]))
    return frame_vectors


def save_charges(charge_positions, charges, filename="out_charges.xyz"):
    print(filename)
    file = open(filename, "w")
    file.write("{}\n".format(len(charge_positions)))
    file.write("s                      x[A]                      y[A]                      z[A]                   "
               "   q[e]\n")

    c = 1
    for xyz, q in zip(charge_positions, charges):
        c += 1
        if q < 0:
            letter = "O"
        else:
            letter = "N"
        file.write("{0:} {1:.16f} {2:.16f} {3:.16f} {4:.16f}\n".format(letter, xyz[0],
                                                                       xyz[1], xyz[2], float(q)))
    file.close()


class ARS():
    def __init__(self, xyz_file_name, pcube, frame_file, pcube_2=None, method="bond", atom_charge_match=None):
        self.method = method
        self.c_positions_local = None
        self.c_positions_global = None
        self.atom_positions = None
        self.atom_positions_plus = None

        # Open XYZ file
        self.c_positions, self.c_charges = read_mdcm_xyz(xyz_file_name)
        self.n_charges = len(self.c_charges)

        # Open Cube files
        self.atom_positions, self.atom_names = read_cube_file(pcube)
        self.n_atoms = len(self.atom_names)

        if pcube_2 is not None:
            self.atom_positions_plus, atom_names = read_cube_file(pcube_2)
            self.n_atoms_2 = len(atom_names)

        #  Test for consistency
        # self.test()

        # Get frames
        self.frame = open(frame_file).readlines()
        self.frame_atoms = []
        self.frames = self.frame[1:]
        self.n_frames = len(self.frames)
        for f in self.frames:
            a1, a2, a3 = f.split()
            self.frame_atoms.append([int(a1), int(a2), int(a3)])

        #  Match charges to closest atoms
        if atom_charge_match is None:
            self.charge_atom_associations, self.atom_charge_dict = self.match_charges()
        else:
            self.read_charge_atom_associations(atom_charge_match)
        print(self.charge_atom_associations)
        print(self.match_charges()[0])

        # Calculate local axes and transform charges
        # Calculate the new axes for each frame
        self.atom_positions = np.array(self.atom_positions)

        if pcube_2 is not None:
            self.atom_positions_plus = np.array(self.atom_positions_plus)

        self.frame_vectors = get_local_axis(self.atom_positions, self.frame_atoms, method=self.method)

        if pcube_2 is not None:
            self.frame_vectors_plus = get_local_axis(self.atom_positions_plus, self.frame_atoms, method=self.method)

        self.c_positions_local = self.global_to_local()

        if pcube_2 is not None:
            self.charge_positions_plus = self.local_to_global()

    def save_charge_atom_associations(self, filename="out.acd"):
        f = open(filename+".acd", "w")
        for charge, atom in self.charge_atom_associations:
            f.write(f"{charge} {atom} \n")
        f.close()

    def read_charge_atom_associations(self, path):
        f = open(path).readlines()
        charge_atom_associations = []
        atom_charge_dict = {}
        for i, line in enumerate(f):
            c, a = line.split()
            c_a = [int(c), int(a)]
            charge_atom_associations.append(c_a)
            c = c_a[0]
            a = c_a[1]
            if a not in atom_charge_dict:
                atom_charge_dict[a] = [c]
            else:
                atom_charge_dict[a].append(c)
        self.charge_atom_associations = charge_atom_associations
        self.atom_charge_dict = atom_charge_dict

    def get_angle(self, a, b, c):
        atoms = self.atom_positions
        p = np.array([atoms[a], atoms[b], atoms[c]])
        return angle(*p)

    def get_dih(self, a, b, c, d):
        atoms = self.atom_positions
        p = np.array([atoms[a], atoms[b], atoms[c], atoms[d]])
        return dihedral2(p)

    def get_dih_2(self, a, b, c, d):
        atoms = self.atom_positions_plus
        p = np.array([atoms[a], atoms[b], atoms[c], atoms[d]])
        return dihedral2(p)

    # def align_in_global(self, filename_template=None):
    #     self.rotation, rmsd = Kabsch.align_vectors(self.atom_positions, self.atom_positions_plus)
    #     self.rotation = self.rotation.as_matrix()
    #     tmp_atom_positions = self.rotation.dot(self.atom_positions.T).T
    #     if filename_template is None:
    #         save_xyz(tmp_atom_positions, self.atom_names)
    #     else:
    #         save_xyz(tmp_atom_positions, self.atom_names, filename=filename_template.format("molecule"))
    #     tmp_charge_positions = self.rotation.dot(self.c_positions.T).T
    #     if filename_template is None:
    #         save_charges(tmp_charge_positions, self.c_charges)
    #     else:
    #         save_charges(tmp_charge_positions, self.c_charges, filename=filename_template.format("charges"))
    #
    #     print(rmsd)

    def get_c_positions_local(self):
        return self.c_positions_local

    def plot1(self):
        plot_labels = False
        plot_pos_1 = True
        plot_pos_2 = True
        plot_vectors = False

        fig = plt.figure()
        ax = Axes3D(fig, elev=0, azim=60)
        # Transpose to the right shape for plotting
        a_p = np.array(self.atom_positions).T
        a_p1 = np.array(self.atom_positions_plus).T
        c_p = np.array(self.c_positions).T
        c_p_l = np.array(self.c_positions_local).T
        c_p_g = np.array(self.c_positions_global).T
        # Plotting axes
        if plot_vectors:
            for frame_vec in self.frame_vectors:
                for il, local_vector_i in enumerate(frame_vec):
                    atom_index = self.frame_atoms[0][il] - 1
                    self.plot_axe(il, local_vector_i, atom_index)
        #  Plotting points
        if plot_pos_1:
            ax.plot(a_p[0], a_p[1], a_p[2], c='gray', linestyle='None', marker="o")
            ax.plot(c_p[0], c_p[1], c_p[2], marker="o", c='orange', linestyle='None', alpha=0.8)
            if plot_labels:
                label = ['{:d}'.format(i) for i in range(self.n_charges)]
                for i, pos in enumerate(c_p.T):
                    ax.text(pos[0], pos[1], pos[2], label[i])
        # ax.plot(c_p_l[0], c_p_l[1], c_p_l[2], marker="o", c="g", linestyle = 'None', alpha=0.8)
        # label = ['{:d}'.format(i) for i in range(n_charges)]
        # for i, pos in enumerate(c_p_l.T):
        #     ax.text(pos[0], pos[1], pos[2], label[i])
        if plot_pos_2:
            ax.plot(a_p1[0], a_p1[1], a_p1[2], c='k', linestyle='None', marker="o")
            ax.plot(c_p_g[0], c_p_g[1], c_p_g[2], marker="x", c="r", linestyle='None', alpha=0.8)
            if plot_labels:
                label = ['{:d}'.format(i) for i in range(self.n_charges)]
                for i, pos in enumerate(c_p_g.T):
                    ax.text(pos[0], pos[1], pos[2], label[i])

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)

        if plot_labels:
            plt.legend()

        plt.show()

    def plot_axe(self, il, local_vector, atom_index, c="k"):
        atom_pos = self.atom_positions[atom_index]
        x = [atom_pos[0], atom_pos[0] + local_vector[0][0]]
        y = [atom_pos[1], atom_pos[1] + local_vector[0][1]]
        z = [atom_pos[2], atom_pos[2] + local_vector[0][2]]
        plt.plot(x, y, z, c='r', label="x")
        x = [atom_pos[0], atom_pos[0] + local_vector[1][0]]
        y = [atom_pos[1], atom_pos[1] + local_vector[1][1]]
        z = [atom_pos[2], atom_pos[2] + local_vector[1][2]]
        plt.plot(x, y, z, '--g', label="y")
        x = [atom_pos[0], atom_pos[0] + local_vector[2][0]]
        y = [atom_pos[1], atom_pos[1] + local_vector[2][1]]
        z = [atom_pos[2], atom_pos[2] + local_vector[2][2]]
        plt.plot(x, y, z, ':b', label="z")
        print("check for orthogality: ", np.dot(local_vector[2], local_vector[0]))

    def match_charges(self):
        # Match each charge to a nucleus
        charge_atom_associations = []
        atom_charge_dict = {}
        for i_charge in range(self.n_charges):
            #  initial distance, which can be compared to find smaller values
            min_distance = np.Inf
            atom_association = None
            for j_atom in range(self.n_atoms):
                d = distance.euclidean(self.c_positions[i_charge], self.atom_positions[j_atom])
                if d < min_distance:
                    atom_association = j_atom
                    min_distance = d
            charge_atom_associations.append([i_charge, atom_association])

            if atom_association not in list(atom_charge_dict.keys()):
                atom_charge_dict[atom_association] = [i_charge]
            else:
                atom_charge_dict[atom_association].append(i_charge)
        return charge_atom_associations, atom_charge_dict

    def global_to_local(self):
        #  Find the position of the charges in the local axes
        #  Create a new array for the 'local' charges
        c_pos_shape = np.array(self.c_positions).shape
        c_positions_local = np.zeros(c_pos_shape)

        used_atoms = []
        for f in range(self.n_frames):
            #  Loop through the atoms in the frame
            for ai, atom_index in enumerate(self.frame_atoms[f]):
                atom_index -= 1
                if atom_index in list(self.atom_charge_dict.keys()) and atom_index not in used_atoms:
                    charges = self.atom_charge_dict[atom_index]
                    ex, ey, ez = self.frame_vectors[f][ai]
                    #  Find the associated charges for that atom, and loop
                    for charge in charges:
                        c_pos_global = self.c_positions[charge]
                        atom_pos_xyz = self.atom_positions[atom_index]

                        #  Find the distance between the charge and the atom it belongs to
                        r = np.array(c_pos_global) - np.array(atom_pos_xyz)

                        local_x_pos = np.dot(ex, r)
                        local_y_pos = np.dot(ey, r)
                        local_z_pos = np.dot(ez, r)

                        c_positions_local[charge][0] = local_x_pos
                        c_positions_local[charge][1] = local_y_pos
                        c_positions_local[charge][2] = local_z_pos

                used_atoms.append(atom_index)

        return c_positions_local

    def save_charges_local(self, output_filename):
        # output_filename_split = output_filename.split("/")
        # output_filename = "local_" + output_filename_split[-1]
        # output_filename = os.path.join(*output_filename_split[:-1], output_filename)
        save_charges(self.c_positions_local,
                     self.c_charges, filename=output_filename + ".local")

    def set_charge_positions_plus(self, charge_positions):
        self.charge_positions_plus = charge_positions

    def set_local_charge_positions(self, charge_positions):
        self.c_positions_local = charge_positions

    def save_charges_global(self, output_filename):
        save_charges(self.charge_positions_plus,
                     self.c_charges, filename=output_filename + ".global")

    def get_distance_charges(self):
        return Kabsch.align_vectors(self.charge_positions_plus, self.c_positions)[1]

    def get_distance_atoms(self):
        return Kabsch.align_vectors(self.atom_positions, self.atom_positions_plus)[1]

    def local_to_global(self):
        #  Find the position of the charges in the local axes
        #  Create a new array for the 'local' charges
        c_pos_shape = np.array(self.c_positions).shape
        c_new_local = np.zeros(c_pos_shape)
        c_positions_global = np.zeros(c_pos_shape)

        used_atoms = []
        for f in range(self.n_frames):
            #  Loop through the atoms in the frame
            for ai, atom_index in enumerate(self.frame_atoms[f]):
                atom_index -= 1
                if atom_index in list(self.atom_charge_dict.keys()) and atom_index not in used_atoms:
                    charges = self.atom_charge_dict[atom_index]
                    ex, ey, ez = self.frame_vectors_plus[f][ai]
                    # Find the associated charges for that atom, and loop
                    for charge in charges:
                        c_pos_local = self.c_positions_local[charge]
                        atom_pos_xyz = self.atom_positions_plus[atom_index]

                        c_l_x = c_pos_local[0]
                        c_l_y = c_pos_local[1]
                        c_l_z = c_pos_local[2]

                        x_vec = np.multiply(ex, c_l_x)
                        y_vec = np.multiply(ey, c_l_y)
                        z_vec = np.multiply(ez, c_l_z)

                        sum_of_components = x_vec + y_vec + z_vec
                        #  translate back to the center of atoms (for the new conformation)
                        c_positions_global[charge] = sum_of_components + atom_pos_xyz

                used_atoms.append(atom_index)
        return c_positions_global


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ARS')
    parser.add_argument('-charges', help='.')
    parser.add_argument('-pcube', help='.')
    parser.add_argument('-pcube2', help='.', default=None)
    parser.add_argument('-frames', help='.')
    parser.add_argument('-output', help='.')
    parser.add_argument('-acd', help='.', default=None)

    args = parser.parse_args()

    ARS_obj = ARS(args.charges, args.pcube, args.frames, pcube_2=args.pcube2, method="bond", atom_charge_match=args.acd)
    ARS_obj.save_charges_global(args.output)
    ARS_obj.save_charges_local(args.output)
    ARS_obj.save_charge_atom_associations(filename=args.output)

    print(f"RMSD_ATOMS = {ARS_obj.get_distance_atoms()}")
    print(f"RMSD_CHARGES = {ARS_obj.get_distance_charges()}")

    # dih = False
    # if len(sys.argv) > 6:
    #     dih = [int(x) for x in sys.argv[6].split("_")]

    # if dih:
    #     dihedral = ARS_obj.get_dih_2(*dih)
    #     print(f"Dihedral {dih} = {dihedral}")
