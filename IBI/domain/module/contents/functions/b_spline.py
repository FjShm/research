import numpy as np
from typing import Union
from module.contents.parameters import PRM
from tqdm import tqdm


# b-spline basis function
def b_basis(
    t: Union[int, float], t_list: list, j: int, k: int
) -> Union[float, np.array]:
    tj = t_list[j]
    tjp1 = t_list[j + 1]
    if k == 0:
        return 1 if tj <= t < tjp1 else 0
    tjk = t_list[j + k]
    tjkp1 = t_list[j + k + 1]
    term1 = (t - tj) / (tjk - tj) * b_basis(t, t_list, j, k - 1)
    term2 = (tjkp1 - t) / (tjkp1 - tjp1) * b_basis(t, t_list, j + 1, k - 1)
    return term1 + term2


# b-spline basis function(differential)
def db_basis(
    t: Union[int, float], t_list: list, j: int, k: int
) -> Union[float, np.array]:
    if k == 0:
        return 0
    tj = t_list[j]
    tjp1 = t_list[j + 1]
    tjk = t_list[j + k]
    tjkp1 = t_list[j + k + 1]
    term1 = b_basis(t, t_list, j, k - 1) / (tjk - tj)
    term2 = (t - tj) / (tjk - tj) * db_basis(t, t_list, j, k - 1)
    term3 = -b_basis(t, t_list, j + 1, k - 1) / (tjkp1 - tjp1)
    term4 = (tjkp1 - t) / (tjkp1 - tjp1) * db_basis(t, t_list, j + 1, k - 1)
    return term1 + term2 + term3 + term4


# b-spline main
def b_spline(
    px: Union[np.array, list], py: Union[np.array, list], n: int = 2, num=0
) -> tuple:
    if len(px) != len(py):
        print("The lengths of the x,y arrays are different.")
        print("(module.contents.functions.b_spline.py)")
        exit()
    lenp = len(px) if num == 0 else num
    m = lenp + n + 1
    t_list = list(np.linspace(0, m - 1, num=m, dtype=int))
    t = np.linspace(t_list[n], t_list[m - n - 1], num=lenp)
    # t = t_list[n:m-n]
    Sx = np.zeros(lenp)
    Sy = np.zeros(lenp)
    dSx = np.zeros(lenp)
    dSy = np.zeros(lenp)
    for ti, t_ in enumerate(tqdm(t)):
        for i in range(m - n - 1):
            if (t_ < t_list[i]) or (t_list[i + n + 1] < t_):
                continue
            bb = b_basis(t_, t_list, i, n)
            dbb = db_basis(t_, t_list, i, n)
            Sx[ti] += px[i] * bb
            Sy[ti] += py[i] * bb
            dSx[ti] += px[i] * dbb
            dSy[ti] += py[i] * dbb
    return list(Sx), list(Sy), dSy / dSx


# b-spline main(differential)
def db_spline(
    px: Union[np.array, list], py: Union[np.array, list], n: int = 2, num: int = 0
) -> list:
    if len(px) != len(py):
        print("The lengths of the x,y arrays are different.")
        print("(module.contents.functions.b_spline.py)")
        exit()
    lenp = len(px) if num == 0 else num
    m = lenp + n + 1
    t_list = list(np.linspace(0, m - 1, num=m, dtype=int))
    t = np.linspace(t_list[n], t_list[m - n - 1], num=lenp)
    Sx = np.zeros(lenp)
    Sy = np.zeros(lenp)
    for ti, t_ in enumerate(tqdm(t)):
        for i in range(m - n - 1):
            if (t_ < t_list[i]) or (t_list[i + n + 1] < t_):
                continue
            Sx[ti] += px[i] * db_basis(t_, t_list, i, n)
            Sy[ti] += py[i] * db_basis(t_, t_list, i, n)
    return Sy / Sx
