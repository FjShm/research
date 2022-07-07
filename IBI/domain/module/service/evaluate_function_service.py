from typing import Union
from module.adapter.evaluate_function_adapter import EvaluateFunctionAdapter
from module.service.get_distribution_service import GetDistribution


class EvaluateFunction:

    __InputData: dict
    dist: GetDistribution
    f: Union[float, None]

    def __init__(self, InputData: dict, dist: GetDistribution) -> None:
        self.__InputData = InputData
        self.dist = dist

    def evaluate(self):
        self.f = EvaluateFunctionAdapter(
            self.dist, self.__InputData["eval_type"], self.__InputData["eval_sigma"]
        ).request()
        return self
