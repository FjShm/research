from typing import Union
from module.facade.get_distribution_facade import GetDistributionFacade
from module.facade.get_distribution_none_facade import GetDistributionNoneFacade


class GetDistributionAdapter:

    path: str
    ratio: float

    def __init__(self, path, ratio):
        self.path = path
        self.ratio = ratio

    def request(self):
        if self.path is None:
            x, P = GetDistributionNoneFacade()
        else:
            x, P = GetDistributionFacade(self)
        return x, P
