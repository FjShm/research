import numpy as np


def nharmonic5(phi, a1, a2, a3, a4, a5):
    """
    Parameters:
    ----------
        phi : float or numpy.array or pandas.DataFrame, [degree]
        a1, a2, ..., a5 : float, [kcal/mol]

    Returns:
    ----------
        U : float or numpy.array or pandas.DataFrame, [kcal/mol]
    """
    Sum = 0
    a = (a1, a2, a3, a4, a5)
    for i in range(len(a)):
        Sum += a[i] * np.cos(np.deg2rad(phi)) ** i
    return Sum


def dnharmonic5(phi, a1, a2, a3, a4, a5):
    Sum = 0
    a = (a1, a2, a3, a4, a5)
    for i in range(len(a)):
        Sum += -i * a[i] * np.sin(np.deg2rad(phi)) ** (i - 1)
    return Sum
