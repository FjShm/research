from typing import Union
from module.adapter.get_distribution_adapter import GetDistributionAdapter


class GetDistribution:

    x_target: list
    x_CG: Union[list, None]
    target: list
    CG: Union[list, None]
    __InputData: dict

    def __init__(self, InputData: dict) -> None:
        self.x_target = []
        self.x_CG = []
        self.target = []
        self.CG = []
        self.__InputData = InputData

    def get_distribution(self):
        self.x_CG, self.CG = GetDistributionAdapter(
            path=self.__InputData["P_CG_path"], ratio=self.__InputData["ratio"],
        ).request()
        self.x_target, self.target = GetDistributionAdapter(
            path=self.__InputData["P_target_path"], ratio=self.__InputData["ratio"],
        ).request()
        return self
