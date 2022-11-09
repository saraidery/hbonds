import numpy as np
import os

from hbonds.hbondanalysis import ClusterHbondAnalyst
from hbonds.hbondanalysis import PBCHbondAnalyst


file_path = os.path.dirname("~/repo/hbonds/hbonds/tests/")
file_in = os.path.join(file_path, "water_cluster.xyz")
print(file_in)
file_out = "water_cluster.out"
print(file_out)
#file_reference = os.path.join(file_path, "water_PBC_reference.out")

HB = ClusterHbondAnalyst(file_in, file_out)
HB.determine_Hbonds()
HB.print_Hbond_for_O(1)
#
#result = filecmp.cmp(file_out, file_reference, shallow=False)
#assert (result == 0)

#directory = "/Users/sarai/repo/Ammonium-Ammonia/calculations/H2O_pimd/xyz/"
#
#for file in os.listdir(directory):
#    filename = os.fsdecode(file)
#    if filename.endswith(".xyz"):
#        print(filename)
#        a = ClusterHbondAnalyst(os.path.join(directory, filename), "water_analysis_wernet_pimd_2.out")
#        a.determine_Hbonds(hbond_condition="wernet2005")
#        a.print_Hbond_for_O(1)
#
