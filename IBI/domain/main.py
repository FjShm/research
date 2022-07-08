import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import yaml
import pprint

from module.service.get_distribution_service import GetDistribution
from module.service.get_Ui_service import GetUi
from module.service.evaluate_function_service import EvaluateFunction
from module.service.calc_Uip1_service import CalcUip1


def argument_parser():
    parser = argparse.ArgumentParser(description="executes BI or IBI.")
    parser.add_argument("-in", "--input", help="input file path", default="input.yaml")
    return parser.parse_args()


def default_input_data():
    InputData = {
        "type": None,
        "function_type": None,
        "eval_type": "bonded",
        "eval_sigma": 1,
        "P_target_path": None,
        "P_CG_path": None,
        "previous_params_filepath": None,
        "min": None,
        "max": None,
        "ratio": 0.01,
        "bounds": (-np.inf, np.inf),
        "output_dir": "step_ip1",
    }
    return InputData


if __name__ == "__main__":
    args = argument_parser()

    # read_input
    with open(args.input, mode="r") as f:
        Input = yaml.safe_load(f)

    print("\n=== Input Info ===\n")
    pprint.pprint(Input)
    for pair_name, InputData_ in Input.items():
        print(f"\n\n === {pair_name} ===")
        InputData = default_input_data()
        for key, val in InputData_.items():
            if val is not None:
                InputData[key] = val
        InputData["section_name"] = f"{InputData['type']}-{pair_name}"

        # get P_target, P_CG
        dist = GetDistribution(InputData).get_distribution()

        # evaluate function
        evaluater = EvaluateFunction(InputData, dist).evaluate()
        if evaluater.f is not None:
            print(f"\nevaluation function type: {InputData['eval_type']}")
            print(f"The evaluation function f for P_i and P_target is {evaluater.f}")

        # get_Ui
        Ui = GetUi(InputData).get_Ui()

        # IBI
        Uip1 = CalcUip1(InputData, dist, Ui).calc_Uip1()

        # save parameters table
        output_dir = InputData["output_dir"]
        os.makedirs(output_dir, exist_ok=True)
        if InputData["function_type"] != "table":
            if InputData["previous_params_filepath"] is not None:
                output_coeff_fname = os.path.basename(
                    InputData["previous_params_filepath"]
                )
            else:
                output_coeff_fname = f"{InputData['section_name']}.param"
            with open(os.path.join(output_dir, output_coeff_fname), mode="w") as f:
                out_str = ""
                for coeff in Uip1.coeff:
                    out_str += f"{coeff} "
                out_str = out_str[:-1]
                f.write(out_str)
        with open(os.path.join(output_dir, f"{pair_name}.table"), mode="w") as f:
            f.write(Uip1.table)

        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111, ylim=(-1.5, 2))
        if InputData["function_type"] != "table":
            if Ui.U is not None:
                ax.plot(Ui.x, Ui.U, color="k", label=r"$U_i$")
            ax.plot(
                Uip1.x_new, Uip1.Uip1, color="g", label=r"$U_{i+1}$", linestyle="dashed"
            )
            ax.plot(
                Uip1.x_new,
                Uip1.Uip1_fitting,
                color="r",
                label=r"$U_{i+1}{\rm -fitting}$",
            )
        else:
            if Ui.U is not None:
                ax.plot(Ui.x, Ui.U, color="k", label=r"$U_i$")
            ax.plot(Uip1.x_new, Uip1.Uip1, color="r", label=r"$U_{i+1}$")
        ax.legend()
        fig.savefig(
            os.path.join(output_dir, f"{InputData['type']}-{pair_name}.png"),
            bbox_inches="tight",
        )

    print("\n\n--- Succeeded!! ---\n")
