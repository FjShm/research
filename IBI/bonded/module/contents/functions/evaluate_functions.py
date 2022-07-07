import numpy as np
from scipy import integrate


def non_bonded(r: list, gr_CG: list, gr_target: list, sigma: float = 1) -> float:
    r = np.array(r)
    gr_CG = np.array(gr_CG)
    gr_target = np.array(gr_target)
    dgr = gr_target - gr_CG
    integrand = np.exp(-r / sigma) * dgr * dgr
    S = integrate.simpson(integrand, r)
    return S


def bonded(x: list, P_CG: list, P_target: list) -> float:
    x = np.array(x)
    P_CG = np.array(P_CG)
    P_target = np.array(P_target)
    dP = P_target - P_CG
    integrand = dP * dP
    S = integrate.simpson(integrand, x)
    return S
