import numpy as np
from typing import Union


def harmonic_angle(
    theta: Union[np.array, float], k: float, theta0: float
) -> Union[np.array, float]:
    """
    k [kcal/mol rad^2]
    theta0 [degree]
    theta [degree]
    """
    dtheta = np.deg2rad(theta - theta0)
    return k * dtheta ** 2


def harmonic_bond(
    r: Union[np.array, float], k: float, r0: float
) -> Union[np.array, float]:
    """
    k [kcal/mol Ang^2]
    r0 [Ang]
    r [Ang]
    """
    dr = r - r0
    return k * dr ** 2


def dharmonic_angle(
    theta: Union[np.array, float], k: float, theta0: float
) -> Union[np.array, float]:
    dtheta = np.deg2rad(theta - theta0)
    return 2 * k * dtheta


def dharmonic_bond(
    r: Union[np.array, float], k: float, r0: float
) -> Union[np.array, float]:
    dr = r - r0
    return 2 * k * dr
