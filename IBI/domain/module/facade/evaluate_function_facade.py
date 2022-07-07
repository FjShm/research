import numpy as np
from scipy.interpolate import interp1d
from module.contents.evaluate_functions import EvaluateFunctions


class EvaluateFunctionFacade:
    def __call__(self, adapter) -> float:
        eval_func = EvaluateFunctions()(adapter.eval_type)

        x_left = max(adapter.x_target[0], adapter.x_CG[0])
        x_right = min(adapter.x_target[-1], adapter.x_CG[-1])
        x = np.linspace(x_left, x_right, num=1000)

        target = interp1d(adapter.x_target, adapter.P_target)(x)
        CG = interp1d(adapter.x_CG, adapter.P_CG)(x)

        if adapter.eval_type == "bonded":
            f = eval_func(x, CG, target)
        elif adapter.eval_type == "non_bonded":
            f = eval_func(x, CG, target, adapter.sigma)
        else:
            print("eval_type is an invalid value!")
            exit()
        return f
