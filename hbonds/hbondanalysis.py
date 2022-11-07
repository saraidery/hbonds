import numpy as np
import os
from abc import ABC
from scipy.spatial import distance_matrix
from hbonds.io import FileHandlerXYZ

from math import pi

from hbonds.periodic_table import (
    symbol_to_Z,
    Z_to_atomic_weight,
)

class HbondAnalyst(ABC):

    """H-bond analyst class, classifies the H-bonds in a
       water geometry (abstract class)
    """

    @abstractmethod
    def __get_xyz_distance(self, i, j):
        pass

    @abstractmethod
    def __compute_distances(self):
        pass


    def determine_Hbonds(self, hbond_condition="theoretical"):
        self.hbond_condition = hbond_condition
        self.__determine_all_Hbonds()

    def print_Hbond_for_O(self, i):
       f = open(self.out_file, "a")
       f.write(f"\nWater H-bond characterization for oxygen number {i}:\n")
       f.write("------------------------------\n")
       f.write(f"File: {self.in_file}\n")
       f.write(f"Atoms: {self.n_atoms} (O: {self.n_O}, H: {self.n_H}, other: {self.n_atoms - self.n_O - self.n_H})\n")


       f.write("\nIndex  Character      OH-distances\n")
       f.write("--------------------------------------\n")
       dOH = self.D[self.O_indices[i-1]][self.bonded_H(self.O_indices[i-1])]
       f.write("{:3d}      D{:d}A{:d}      {:.4f}       {:.4f}\n".format(self.O_indices[i-1] + 1,
                                                                      self.donating[i-1],
                                                                      self.accepting[i-1],
                                                                      dOH[0],
                                                                      dOH[1]
                                                                     )
        )

       f.write("--------------------------------------\n")
       f.close()



    def __compute_distances(self):
        """Distances between atoms in xyz-geometry with PBC"""

        # Distances in x,y,z-directions between every atom
        if self.PBC:
            xyz_distances = self.xyz[:,None] - self.xyz[None,:]
            self.D = np.zeros([np.size(self.Z), np.size(self.Z)])

            for i in range(np.size(xyz_distances,0)):

                xyz = xyz_distances[i]

                S = np.matmul(xyz, self.toABC)
                S = S - np.rint(S)
                xyz = np.matmul(S, self.fromABC)

                self.D[i] = np.sqrt(np.sum(xyz**2, axis=1))
        else:
            self.D = distance_matrix(self.xyz, self.xyz)


    def __determine_all_Hbonds_theoretical(self):

        accepting = np.zeros(self.n_atoms, dtype=int)
        donating = np.zeros(self.n_atoms, dtype=int)

        A = self.non_bonded_neighbours

        for i in self.O_indices:

            for j in self.O_indices:

                if (A[i][j] == 0):
                    continue

                xOO, yOO, zOO = self.__get_xyz_distance(i, j)
                dOO = self.D[i][j]

                for k in self.bonded_H(i):
                    xOH, yOH, zOH = self.__get_xyz_distance(i, k)
                    dOH = self.D[i][k]

                    cos_theta = (xOO*xOH + yOO*yOH + zOO*zOH)/(dOO*dOH)
                    theta = radians_to_degrees(np.arccos(cos_theta))

                    if (theta < 30):
                        accepting[j] = accepting[j] + 1
                        donating[i] = donating[i] + 1

        self.accepting = accepting[self.O_indices]
        self.donating = donating[self.O_indices]


    def __determine_all_Hbonds_wernet(self):
        accepting = np.zeros(self.n_atoms, dtype=int)
        donating = np.zeros(self.n_atoms, dtype=int)

        A = self.non_bonded_neighbours

        for i in self.O_indices:

            for j in self.O_indices:

                if (A[i][j] == 0):
                    continue

                xOO, yOO, zOO = self.__get_xyz_distance(i, j)
                dOO = self.D[i][j]

                for k in self.bonded_H(i):
                    xOH, yOH, zOH = self.__get_xyz_distance(i, k)
                    dOH = self.D[i][k]

                    cos_theta = (xOO*xOH + yOO*yOH + zOO*zOH)/(dOO*dOH)
                    theta = radians_to_degrees(np.arccos(cos_theta))

                    if (theta < 45):
                        if (dOO < 3.3 - 0.00044*theta**2):
                            accepting[j] = accepting[j] + 1
                            donating[i] = donating[i] + 1

        self.accepting = accepting[self.O_indices]
        self.donating = donating[self.O_indices]

    def __determine_all_Hbonds(self):

        if self.hbond_condition == "theoretical":
            self.__determine_all_Hbonds_theoretical()
        elif self.hbond_condition == "wernet2005":
            self.__determine_all_Hbonds_wernet()
        else:
            print("Hbond condition should be either theoretical or wernet2005, default is theoretical")

    def print_summary(self):

        f = open(self.out_file, "a")
        f.write("\nWater H-bond characterization:\n")
        f.write("------------------------------\n")
        f.write(f"File: {self.in_file}\n")
        f.write(f"Atoms: {self.n_atoms} (O: {self.n_O}, H: {self.n_H}, other: {self.n_atoms - self.n_O - self.n_H})\n")
        f.write(f"\nH-bond condition: {self.hbond_condition}\n")
        f.write("\nIndex  Character      OH-distances\n")
        f.write("--------------------------------------\n")
        for i in range(self.n_O):
            dOH = self.D[self.O_indices[i]][self.bonded_H(self.O_indices[i])]
            f.write("{:3d}      D{:d}A{:d}      {:.4f}       {:.4f}\n".format(self.O_indices[i] + 1,
                                                                           self.donating[i],
                                                                           self.accepting[i],
                                                                           dOH[0],
                                                                           dOH[1]
                                                                          )
                )

        f.write("--------------------------------------\n")
        f.close()

    @property
    def n_O(self):
        return (self.Z == 8).sum()

    @property
    def O_indices(self):
        return np.where(self.Z == 8)[0]

    @property
    def n_H(self):
        return (self.Z == 1).sum()

    @property
    def H_indices(self):
        return np.where(self.Z == 1)[0]

    @property
    def n_atoms(self):
        return np.size(self.Z)

    @property
    def bond_matrix(self):
        B = np.zeros([self.n_atoms, self.n_atoms], dtype=int)

        B[np.where((self.D < 1.3) & (self.D > 0.0))] = 1
        return B

    @property
    def non_bonded_neighbours(self):

        B = np.zeros([self.n_atoms, self.n_atoms], dtype=int)
        B[np.where((self.D < 3.5) & (self.D > 1.3))] = 1
        return B

    @property
    def bonds_per_atom(self):
        return np.sum(self.bond_matrix, axis=1)

    @property
    def bonds_per_O(self):
        return self.bonds_per_atom[self.O_indices]

    @property
    def bonds_per_H(self):
        return self.bonds_per_atom[self.H_indices]

    @property
    def bonded_H(self, i):
        return np.where(self.bond_matrix[i] == 1)[0]


class PBCHbondAnalyst(HbondAnalyst):

    """H-bond analyst class, classifies the H-bonds in a
       water geometry (periodic box)
    """

    def __init__(self, in_file, out_file, a, b, c, alpha=90, beta=90, gamma=90, angle_unit="degrees"):

        """
        """
        self.in_file = in_file
        self.out_file = out_file

        xyz_reader = FileHandlerXYZ(self.in_file)
        self.symbols, self.xyz = xyz_reader.read()
        self.Z = np.fromiter(map(symbol_to_Z, np.atleast_1d(self.symbols)), dtype=int)

        self.__set_PBC(a, b, c, alpha, beta, gamma, angle_unit)
        self.__compute_distances()

    def __set_PBC(self, a, b, c, alpha=90, beta=90, gamma=90, angle_unit="degrees"):

        self.PBC = True

        # Cell parameters
        self.a = a
        self.b = b
        self.c = c

        if (angle_unit == "degrees"):
            self.alpha = degrees_to_radians(alpha)
            self.beta  = degrees_to_radians(beta)
            self.gamma = degrees_to_radians(gamma)


    def __fromABC(self):

        cos_a = np.cos(self.alpha)
        cos_b = np.cos(self.beta)
        cos_c = np.cos(self.gamma)
        sin_c = np.sin(self.gamma)

        T = np.zeros([3,3])
        T[0][0] = self.a
        T[1][0] = self.b*cos_c
        T[1][1] = self.b*sin_c
        T[2][0] = self.c*cos_b
        T[2][1] = self.c*(cos_a - cos_b*cos_c)/sin_c
        T[2][2] = self.c*np.sqrt((1.0 + 2.0*cos_a*cos_b*cos_c - cos_a**2 - cos_b**2 - cos_c**2)/sin_c)

        return T

    def __toABC(self):

        U = np.zeros([3,3])
        T = self.fromABC
        U[0][0] = 1.0/T[0][0]
        U[1][0] = -T[1][0]/(T[0][0]*T[1][1])
        U[1][1] = 1.0/T[1][1]
        U[2][0] = (-T[1][1]*T[2][0]+T[1][0]*T[2][1])/(T[0][0]*T[1][1]*T[2][2])
        U[2][1] = -T[2][1]/(T[1][1]*T[2][2])
        U[2][2] = 1.0/T[2][2]

        return U

    def __get_xyz_distance(self, i, j):

        xyz =  self.xyz[i,:] - self.xyz[j,:]

        ABC = np.matmul(xyz, self.toABC)
        ABC = ABC - np.rint(ABC)
        xyz = np.matmul(ABC, self.fromABC)

        return xyz[0], xyz[1], xyz[2]


    def __compute_distances(self):
        """Distances between atoms in xyz-geometry with PBC"""

        # Distances in x,y,z-directions between every atom
        xyz_distances = self.xyz[:,None] - self.xyz[None,:]
        self.D = np.zeros([np.size(self.Z), np.size(self.Z)])

        for i in range(np.size(xyz_distances,0)):

            xyz = xyz_distances[i]

            S = np.matmul(xyz, self.toABC)
            S = S - np.rint(S)
            xyz = np.matmul(S, self.fromABC)

            self.D[i] = np.sqrt(np.sum(xyz**2, axis=1))


class ClusterHbondAnalyst(HbondAnalyst):

    """H-bond analyst class, classifies the H-bonds in a
       water geometry (periodic box)
    """

    def __init__(self, in_file, out_file):

        """
        """
        self.in_file = in_file
        self.out_file = out_file

        xyz_reader = FileHandlerXYZ(self.in_file)
        self.symbols, self.xyz = xyz_reader.read()
        self.Z = np.fromiter(map(symbol_to_Z, np.atleast_1d(self.symbols)), dtype=int)

        self.__compute_distances()

    def __get_xyz_distance(self, i, j):

        xyz =  self.xyz[i,:] - self.xyz[j,:]

        return xyz[0], xyz[1], xyz[2]

    def __compute_distances(self):
        """Distances between atoms in xyz-geometry"""

        self.D = distance_matrix(self.xyz, self.xyz)


def degrees_to_radians(angle):
    return angle*pi/180.0

def radians_to_degrees(angle):
    return angle*180.0/pi


