import numpy as np
from scipy.interpolate import splrep, splev
from typing import Union
from module.contents.parameters import PRM


# b-spline
def b_spline_scipy(
    px: Union[np.array, list], py: Union[np.array, list], num: int = 100, spacing=False,
) -> list:
    if len(px) != len(py):
        print("The lengths of the x,y arrays are different.")
        print("(module.contents.functions.b_spline.py)")
        exit()

    # fitting/smoothing
    n_interior_knots = num
    qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
    knots = np.quantile(px, qs)
    tck = splrep(px, py, t=knots, k=3)
    if spacing:
        px = np.linspace(px[0], px[-1], num=len(px))
        py_smooth = list(splev(px, tck))
        return px, py_smooth
    else:
        py_smooth = list(splev(px, tck))
    return py_smooth


# b-spline main(differential)
def db_spline_scipy():
    pass
