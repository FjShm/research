from module.contents.functions.harmonic import (
    harmonic_bond,
    harmonic_angle,
    dharmonic_bond,
    dharmonic_angle,
)
from module.contents.functions.gaussian import gaussian, dgaussian
from module.contents.functions.nharmonic5 import nharmonic5, dnharmonic5
from module.contents.functions.nharmonic7 import nharmonic7, dnharmonic7


class PotentialFunctions:
    def __call__(self, func_name: str, d=False):
        if func_name == "harmonic_angle":
            if d:
                return dharmonic_angle
            else:
                return harmonic_angle
        if func_name == "harmonic_bond":
            if d:
                return dharmonic_bond
            else:
                return harmonic_bond
        elif func_name == "gaussian":
            if d:
                return dgaussian
            else:
                return gaussian
        elif func_name == "nharmonic5":
            if d:
                return dnharmonic5
            else:
                return nharmonic5
        elif func_name == "nharmonic7":
            if d:
                return dnharmonic7
            else:
                return nharmonic7
        else:
            print(f"The given potential function type {func_name} is invalid.")
            exit()
