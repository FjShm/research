import numpy as np
from module.contents.parameters import PRM


def gaussian(theta, a1, w1, theta1, a2, w2, theta2) -> float:
    """
                              a1                   theta - theta1                 a2                   theta - theta2
    U(theta) = -kBT ln[ --------------- exp{ -2 ( ---------------- )**2 } + --------------- exp{ -2 ( ---------------- )**2 } ]
                         w1*sqrt(pi/2)                   w1                  w2*sqrt(pi/2)                   w2

    units for LAMMPS units "real"

    Parameters:
    ----------
        theta : float or numpy.array or pandas.DataFrame, [degree]
        a1,a2,w1,w2,theta1,theta2 : float

    Returns:
    ----------
        U : float or numpy.array or pandas.DataFrame, [kcal/mol]
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


def dgaussian(theta, a1, w1, theta1, a2, w2, theta2) -> float:
    return theta * 0
