from typing import Union
from module.adapter.calc_Uip1_adapter import CalcUip1Adapter
from module.service.get_distribution_service import GetDistribution
from module.service.get_Ui_service import GetUi


class CalcUip1:

    x_new: list
    Uip1: list
    Uip1_fitting: Union[list, None]
    coeff: list
    table: str
    __InputData: dict

    def __init__(self, InputData: dict, dist: GetDistribution, Ui: GetUi) -> None:
        self.__dist = dist
        self.__Ui = Ui
        self.__InputData = InputData

    def calc_Uip1(self):
        (
            self.x_new,
            self.Uip1,
            self.Uip1_fitting,
            self.coeff,
            self.table,
        ) = CalcUip1Adapter(
            function_type=self.__InputData["function_type"],
            section_name=self.__InputData["section_name"],
            bounds=tuple(self.__InputData["bounds"]),
            dist=self.__dist,
            Ui=self.__Ui,
            Min=self.__InputData["min"],
            Max=self.__InputData["max"],
            num=self.__InputData["num_table"],
            ibi_accelerator=self.__InputData["ibi_accelerator"],
            shift=self.__InputData["shift_U_min_to_zero"],
        ).request()
        return self
