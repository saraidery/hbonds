import numpy as np
import os

from hbonds.hbondanalysis import ClusterHbondAnalyst

directory = "/Users/sarai/repo/Ammonium-Ammonia/calculations/H2O_pimd/xyz/"

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".xyz"):
        print(filename)
        a = ClusterHbondAnalyst(os.path.join(directory, filename), "water_analysis_wernet_pimd_2.out")
        a.determine_Hbonds(hbond_condition="wernet2005")
        a.print_Hbond_for_O(1)

