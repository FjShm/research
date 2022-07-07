import numpy as np
from scipy.interpolate import interp1d

from typing import Union
from module.facade.calc_Uip1_table_facade import CalcUip1TableFacade
from module.facade.calc_Uip1_function_facade import CalcUip1FunctionFacade


class CalcUip1Adapter:

    P_target: list
    P_CG: list
    Ui: list

    def __init__(
        self,
        function_type: str,
        section_name: str,
        bounds: tuple,
        x_target: list,
        x_CG: Union[list, None],
        x_Ui: Union[list, None],
        P_target: list,
        P_CG: Union[list, None],
        Ui: Union[list, None],
        Min: Union[float, None],
        Max: Union[float, None],
    ) -> None:
        self.function_type = function_type
        self.section_name = section_name
        self.P_target = P_target
        self.P_CG = P_CG
        self.Ui = Ui
        self.x_target = x_target
        self.x_CG = x_CG
        self.x_Ui = x_Ui
        self.bounds = bounds
        self.Min = Min
        self.Max = Max

    def request(self):
        # x-axisが揃うように線形補間
        if self.P_CG is None:
            self.x_new = self.x_target
        else:
            range_left = max(self.x_target[0], self.x_CG[0], self.x_Ui[0])
            range_right = min(self.x_target[-1], self.x_CG[-1], self.x_Ui[-1])
            num = max(len(self.x_target), len(self.x_CG), len(self.x_Ui))

            # linear interpolation
            lin_target = interp1d(self.x_target, self.P_target)
            lin_CG = interp1d(self.x_CG, self.P_CG)
            lin_Ui = interp1d(self.x_Ui, self.Ui)

            x = np.linspace(range_left, range_right, num=num)
            self.P_target = list(lin_target(x))
            self.P_CG = list(lin_CG(x))
            self.Ui = list(lin_Ui(x))
            self.x_new = x

        if self.function_type == "table":
            Uip1, Uip1_fitting, coeff, table = CalcUip1TableFacade()(self)
        else:
            Uip1, Uip1_fitting, coeff, table = CalcUip1FunctionFacade()(self)
        return self.x_new, Uip1, Uip1_fitting, coeff, table
