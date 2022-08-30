import numpy as np
from typing import Union
from module.contents.parameters import PRM


class Flatten:

    a = PRM.mypoly_a
    x0 = PRM.mypoly_x0

    def tanh(self, x):
        return 0.5 - np.tanh(self.a * (x - self.x0)) * 0.5

    def dtanh(self, x) -> float:
        cosh = np.cosh(self.a * (x - self.x0))
        return -0.5 * self.a / (cosh * cosh)


def pair_myPolynominal(r: Union[np.array, float], *args) -> Union[np.array, float]:
    rcut = PRM.rcut
    fUip1 = np.poly1d(args)
    Uip1 = fUip1(r - rcut) * Flatten().tanh(r)
    return Uip1


def dpair_myPolynominal(r: Union[np.array, float], *args) -> Union[np.array, float]:
    rcut = PRM.rcut
    tanh = Flatten().tanh(r)
    dtanh = Flatten().dtanh(r)
    fUip1 = np.poly1d(args)
    fdUip1 = np.polyder(fUip1)
    Uip1 = fUip1(r - rcut)
    dUip1 = fdUip1(r - rcut) * tanh + Uip1 * dtanh
    return dUip1
