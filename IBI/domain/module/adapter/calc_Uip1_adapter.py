import numpy as np
from scipy.interpolate import interp1d

from typing import Union
from module.facade.calc_Uip1_table_facade import CalcUip1TableFacade
from module.facade.calc_Uip1_polynominal_facade import CalcUip1PolynominalFacade
from module.facade.calc_Uip1_myPolynominal_facade import CalcUip1MyPolynominalFacade
from module.facade.calc_Uip1_function_facade import CalcUip1FunctionFacade
from module.service.get_distribution_service import GetDistribution
from module.service.get_Ui_service import GetUi


class CalcUip1Adapter:

    P_target: list
    P_CG: Union[list, None]
    Ui: Union[list, None]
    x_target: list
    x_CG: Union[list, None]
    x_Ui: Union[list, None]

    def __init__(
        self,
        function_type: str,
        section_name: str,
        bounds: tuple,
        dist: GetDistribution,
        Ui: GetUi,
        Min: float,
        Max: float,
        num: int,
        ibi_accelerator: float = 1,
    ) -> None:
        self.function_type = function_type
        self.section_name = section_name
        self.P_target = dist.target
        self.P_CG = dist.CG
        self.Ui = Ui.U
        self.x_target = dist.x_target
        self.x_CG = dist.x_CG
        self.x_Ui = Ui.x
        self.bounds = bounds
        self.Min = Min
        self.Max = Max
        self.num = num
        self.ibi_accelerator = ibi_accelerator

    def request(self):
        # x-axisが揃うように線形補間
        if self.P_CG is None:
            if self.num == 0:
                self.x_new = self.x_target
            else:
                num = self.num
                # linear interpolation
                lin_target = interp1d(self.x_target, self.P_target)
                # lin_Ui = interp1d(self.x_Ui, self.Ui)

                x = np.linspace(self.x_target[0], self.x_target[-1], num=num)
                self.P_target = list(lin_target(x))
                # self.Ui = list(lin_Ui(x))
                self.x_new = x
        else:
            range_left = max(self.x_target[0], self.x_CG[0], self.x_Ui[0])
            range_right = min(self.x_target[-1], self.x_CG[-1], self.x_Ui[-1])
            if self.num == 0:
                num = max(len(self.x_target), len(self.x_CG), len(self.x_Ui))
            else:
                num = self.num

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
        elif self.function_type == "pair_polynominal":
            Uip1, Uip1_fitting, coeff, table = CalcUip1PolynominalFacade()(self)
        elif self.function_type == "pair_myPolynominal":
            Uip1, Uip1_fitting, coeff, table = CalcUip1MyPolynominalFacade()(self)
        else:
            Uip1, Uip1_fitting, coeff, table = CalcUip1FunctionFacade()(self)
        return self.x_new, Uip1, Uip1_fitting, coeff, table
