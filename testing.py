import numpy as np
import os

from hbonds import io
from hbonds.hbondanalysis import HbondAnalyst

directory = "/Users/sarai/repo/Ammonium-Ammonia/"

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".xyz"):
        a = HbondAnalyst(os.path.join(directory, filename), "water_analysis_wernet.out")
        a.set_PBC(9.85, 9.85, 9.85)
        a.determine_Hbonds(hbond_condition="wernet2005")
        a.print_summary()
