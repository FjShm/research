import numpy as np
from typing import Union
from module.contents.parameters import PRM


def gaussian(
    theta: Union[float, np.array],
    a1: float,
    w1: float,
    theta1: float,
    a2: float,
    w2: float,
    theta2: float,
) -> Union[float, np.array]:
    """
    Parameters:
    ----------
        theta : [degree]
        a1,a2,w1,w2,theta1,theta2 : [degree] or [rad]

    Returns:
    ----------
        U : [kcal/mol]
    """
    kBT = PRM.kBT
    a = (a1, a2)
    w = (w1, w2)
    dthetai = (np.deg2rad(theta - theta1), np.deg2rad(theta - theta2))
    Sum = 0
    for i in range(2):
        Sum += (
            a[i] / (w[i] * np.sqrt(np.pi * 0.5)) * np.exp(-2 * (dthetai[i] / w[i]) ** 2)
        )
    return -kBT * np.log(Sum)


def dgaussian(
    theta: Union[float, np.array],
    a1: float,
    w1: float,
    theta1: float,
    a2: float,
    w2: float,
    theta2: float,
) -> Union[float, np.array]:
    # not yet
    return theta * 0
