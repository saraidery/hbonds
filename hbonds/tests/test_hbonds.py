import numpy as np
import pytest
import os
import filecmp

from hbonds.hbondanalysis import ClusterHbondAnalyst
from hbonds.hbondanalysis import PBCHbondAnalyst

class TestHbonds:
    def test_PBC_geometry(self):

        file_path = os.path.dirname(__file__)
        file_in = os.path.join(file_path, "water_PBC.xyz")
        file_out = os.path.join(file_path, "water_PBC.out")
        file_reference = os.path.join(file_path, "water_PBC_ref.out")

        HB = PBCHbondAnalyst(file_in, file_out, 9.85, 9.85, 9.85)
        HB.determine_Hbonds()
        HB.print_summary()

        ref_list = get_file_list(file_reference)
        out_list = get_file_list(file_out)

        os.remove(file_out)

        assert ref_list == out_list


    def test_PBC_geometry_wernet(self):

        file_path = os.path.dirname(__file__)
        file_in = os.path.join(file_path, "water_PBC.xyz")
        file_out = os.path.join(file_path, "water_PBC.out")
        file_reference = os.path.join(file_path, "water_PBC_ref_wernet.out")

        HB = PBCHbondAnalyst(file_in, file_out, 9.85, 9.85, 9.85)
        HB.determine_Hbonds(hbond_condition="wernet2005")
        HB.print_summary()

        ref_list = get_file_list(file_reference)
        out_list = get_file_list(file_out)

        os.remove(file_out)

        assert ref_list == out_list

    def test_cluster_geometry(self):

        file_path = os.path.dirname(__file__)
        file_in = os.path.join(file_path, "water_cluster.xyz")
        file_out = os.path.join(file_path, "water_cluster.out")
        file_reference = os.path.join(file_path, "water_cluster_ref.out")

        HB = ClusterHbondAnalyst(file_in, file_out)
        HB.determine_Hbonds()
        HB.print_Hbond_for_O(1)

        ref_list = get_file_list(file_reference)
        out_list = get_file_list(file_out)

        os.remove(file_out)

        assert ref_list == out_list


def get_file_list(file_name):
    list_ = []
    with open(file_name, mode="r") as f:
        for line in f.readlines():
            if "File:" not in line:
                list_.append(line)

    return list_
