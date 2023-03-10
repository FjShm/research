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
            Min_parabola_axis=self.__InputData["min_parabola_axis"],
            Max=self.__InputData["max"],
            Max_parabola_axis=self.__InputData["max_parabola_axis"],
            num=self.__InputData["num_table"],
            ibi_accelerator=self.__InputData["ibi_accelerator"],
            Mind_coeff=self.__InputData["mind_coeff"],
            Maxd_coeff=self.__InputData["maxd_coeff"],
            extrapolate_type=self.__InputData["extrapolate_type"],
            shift=self.__InputData["shift_U_min_to_zero"],
        ).request()
        return self
