from typing import Union
from module.adapter.evaluate_function_adapter import EvaluateFunctionAdapter
from module.service.get_distribution_service import GetDistribution


class EvaluateFunction:

    InputData: dict
    dist: GetDistribution
    f: Union[float, None]

    def __init__(self, InputData: dict, dist: GetDistribution):
        self.InputData = InputData
        self.dist = dist

    def evaluate(self):
        self.f = EvaluateFunctionAdapter(
            self.dist, self.InputData["eval_type"], self.InputData["eval_sigma"]
        ).request()
        return self
