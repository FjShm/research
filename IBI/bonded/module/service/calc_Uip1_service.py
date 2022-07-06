from typing import Union
from module.adapter.calc_Uip1_adapter import CalcUip1Adapter
from module.service.get_distribution_service import GetDistribution
from module.service.get_Ui_service import GetUi


class CalcUip1:

    __P_target: list
    __P_CG: list
    __x_target: list
    __x_CG: list
    __Ui: list
    __x_Ui: list

    x_new: list
    Uip1: list
    Uip1_fitting: Union[list, None]
    coeff: list
    table: str
    InputData: dict

    def __init__(self, InputData: dict, dist: GetDistribution, Ui: GetUi) -> None:
        self.__P_target = dist.target
        self.__P_CG = dist.CG
        self.__x_target = dist.x_target
        self.__x_CG = dist.x_CG
        self.__Ui = Ui.U
        self.__x_Ui = Ui.x
        self.InputData = InputData

    def calc_Uip1(self) -> None:
        (
            self.x_new,
            self.Uip1,
            self.Uip1_fitting,
            self.coeff,
            self.table,
        ) = CalcUip1Adapter(
            function_type=self.InputData["function_type"],
            bounds=tuple(self.InputData["bounds"]),
            x_target=self.__x_target,
            x_CG=self.__x_CG,
            x_Ui=self.__x_Ui,
            P_target=self.__P_target,
            P_CG=self.__P_CG,
            Ui=self.__Ui,
        ).request()
