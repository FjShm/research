import numpy as np


def nharmonic7(phi, a1, a2, a3, a4, a5, a6, a7):
    """
    Parameters:
    ----------
        phi : float or numpy.array or pandas.DataFrame, [degree]
        a1, a2, ..., a7 : float, [kcal/mol]

    Returns:
    ----------
        U : float or numpy.array or pandas.DataFrame, [kcal/mol]
    """
    Sum = 0
    a = (a1, a2, a3, a4, a5, a6, a7)
    for i in range(len(a)):
        Sum += a[i] * np.cos(np.deg2rad(phi)) ** i
    return Sum


def dnharmonic7(phi, a1, a2, a3, a4, a5, a6, a7):
    Sum = 0
    a = (a1, a2, a3, a4, a5, a6, a7)
    for i in range(len(a)):
        Sum += -i * a[i] * np.sin(np.deg2rad(phi)) ** (i - 1)
