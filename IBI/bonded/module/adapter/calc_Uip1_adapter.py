import numpy as np

from typing import Union
from module.facade.calc_Uip1_table_facade import CalcUip1TableFacade
from module.facade.calc_Uip1_function_facade import CalcUip1FunctionFacade


class CalcUip1Adapter:

    P_target: list
    P_CG: list
    Ui: list

    def __init__(self, function_type: str, bounds: tuple, x_target: list, x_CG: Union[list, None], x_Ui: Union[list, None], P_target: list, P_CG: Union[list, None], Ui: Union[list, None]):
        self.function_type = function_type
        self.P_target = P_target
        self.P_CG = P_CG
        self.Ui = Ui
        self.x_target = x_target
        self.x_CG = x_CG
        self.Ui = Ui
        self.bounds = bounds

    def request(self):
        # x-axisが揃うように線形補間
        if self.P_CG is None:
            range_left = max(self.x_target[0], self.x_CG[0], self.x_Ui[0])
            range_right = min(self.x_target[-1], self.x_CG[-1], self.x_Ui[-1])
            num = max(len(self.x_target), len(self.x_CG), len(self.x_Ui))

            # linear interpolation
            x = np.linspace(range_left, range_right, num=num)
            lin_target = interp1d(self.x_target, self.P_target)
            lin_CG = interp1d(self.x_CG, self.P_CG)
            lin_Ui = interp1d(self.x_Ui, self.Ui)

            self.P_target = list(lin_target(x))
            self.P_CG = list(lin_CG(x))
            self.Ui = list(lin_Ui(x))
            self.x_new = x_new
        else:
            self.x_new = self.x_target

        if self.function_type == "table":
            Uip1, Uip1_fitting, coeff, table = CalcUip1TableFacade(self)
        else:
            Uip1, Uip1_fitting, coeff, table = CalcUip1FunctionFacade(self)
        return self.x_new, Uip1, Uip1_fitting, coeff, table
