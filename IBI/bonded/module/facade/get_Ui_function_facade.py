import numpy as np
from module.contents.potential_functions import PotentialFunctions


class GetUiFunctionFacade:
    def __call__(self, Ui_adapter):
        function_type = Ui_adapter.function_type
        previous_params_filepath = Ui_adapter.previous_params_filepath
        Min = Ui_adapter.Min
        Max = Ui_adapter.Max
        with open(previous_params_filepath, mode="r") as f:
            for line in f:
                coeff = [float(val) for val in line.split(" ")]
                break
        x = np.linspace(Min, Max, num=500)
        func = PotentialFunctions(function_type)
        U = func(x, *coeff)
        return list(x), list(U)
