import numpy as np
from typing import Union
from module.contents.parameters import PRM


def pair_polynominal(r: Union[np.array, float], *args) -> Union[np.array, float]:
    rcut = PRM.rcut
    fUip1 = np.poly1d(args)
    Uip1 = fUip1(r - rcut)
    return Uip1


def dpair_polynominal(r: Union[np.array, float], *args) -> Union[np.array, float]:
    rcut = PRM.rcut
    fUip1 = np.poly1d(args)
    fdUip1 = np.polyder(fUip1)
    dUip1 = fdUip1(r - rcut)
    return dUip1
