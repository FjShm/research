from module.facade.evaluate_function_facade import EvaluateFunctionFacade
from module.facade.evaluate_function_none_facade import EvaluateFunctionNoneFacade
from module.service.get_distribution_service import GetDistribution


class EvaluateFunctionAdapter:

    dist: GetDistribution
    eval_type: str
    sigma: float

    def __init__(
        self, dist: GetDistribution, eval_type: str = "bonded", sigma: float = 1
    ) -> None:
        self.P_target = dist.target
        self.P_CG = dist.CG
        self.x_target = dist.x_target
        self.x_CG = dist.x_CG
        self.eval_type = eval_type
        self.sigma = sigma

    def request(self) -> float:
        if self.P_CG is None:
            f = EvaluateFunctionNoneFacade()()
        else:
            f = EvaluateFunctionFacade()(self)
        return f
