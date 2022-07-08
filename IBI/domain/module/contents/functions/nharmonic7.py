import numpy as np
from typing import Union


def nharmonic7(
    phi: Union[float, np.array],
    a1: float,
    a2: float,
    a3: float,
    a4: float,
    a5: float,
    a6: float,
    a7: float,
) -> Union[float, np.array]:
    """
    Parameters:
    ----------
        phi : [degree]
        a1, a2, ..., a7 : [kcal/mol]

    Returns:
    ----------
        U : [kcal/mol]
    """
    Sum = 0
    a = (a1, a2, a3, a4, a5, a6, a7)
    for i in range(len(a)):
        Sum += a[i] * np.cos(np.deg2rad(phi)) ** i
    return Sum


def dnharmonic7(
    phi: Union[float, np.array],
    a1: float,
    a2: float,
    a3: float,
    a4: float,
    a5: float,
    a6: float,
    a7: float,
) -> Union[float, np.array]:
    Sum = 0
    a = (a1, a2, a3, a4, a5, a6, a7)
    for i in range(len(a)):
        Sum += -i * a[i] * np.sin(np.deg2rad(phi)) ** (i - 1)
