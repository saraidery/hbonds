import numpy as np
import os


from hbonds import io
from hbonds.hbondanalysis import HbondAnalyst


a = HbondAnalyst("hbonds/tests/H2O_pimd_S01.xyz", 9.85, 9.85, 9.85)
print(a.distances)
print(a.bond_matrix)
print(a.close_O)
