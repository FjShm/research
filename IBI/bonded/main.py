import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import yaml
import pprint

from module.service.get_distribution_service import GetDistribution
from module.service.get_Ui_service import GetUi
from module.service.calc_Uip1_service import CalcUip1

from module.contents.parameters import PRM

from module.common.LAMMPS_potential_table_IO import LAMMPSPotentialTableIO


def argument_parser():
    parser = argparse.ArgumentParser(description="executes BI or IBI.")
    parser.add_argument("-in", "--input", help="input file path", default="input.yaml")
    return parser.parse_args()


def default_input_data():
    InputData = {
        "type": None,
        "function_type": None,
        "P_target_path": None,
        "P_CG_path": None,
        "previous_params_filepath": None,
        "min": None,
        "max": None,
        "ratio": 0.01,
    }
    return InputData


if __name__ == "__main__":
    args = argument_parser()

    # read_input
    with open(args.input, mode="r") as f:
        Input = yaml.safe_load(f)

    for pair_name, InputData_ in Input.items():
        InputData = default_input_data()
        for key, val in InputData_.items():
            if val is not None:
                InputData[key] = val
        pprint.pprint(InputData)
        exit()

        # get P_target, P_CG
        dist = GetDistribution(InputData).get_distribution()

        # get_Ui
        Ui = GetUi(InputData).get_Ui()

        # IBI
        Uip1 = CalcUip1(InputData, dist, Ui).calc_Uip1()

        # save parameters table

        # plot
