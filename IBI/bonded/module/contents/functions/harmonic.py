import numpy as np


def harmonic(theta, k: float, theta0: float) -> float:
    """
    k [kcal/mol rad^2]
    theta0 [degree]
    theta [ degree ]
    """
    dtheta = np.deg2rad(theta - theta0)
    return k * dtheta ** 2


def dharmonic(theta, k, theta0) -> float:
    dtheta = np.deg2rad(theta - theta0)
    return 2 * k * dtheta
