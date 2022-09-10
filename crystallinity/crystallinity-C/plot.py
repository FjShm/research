import sys
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def func(x, a: float, b: float, idx: float):
    return a * np.abs(x)**idx + b

class PRM:
    fs2ns = 1e-6


if __name__ == "__main__":
    txt = np.loadtxt(sys.argv[1], dtype="str", delimiter=" ")
    erate = float(input("please type erate [ns^-1]: "))
    dt = float(input("type dt [fs]: "))
    fig = plt.figure(dpi=600)
    ax = fig.add_subplot(
        111, xlabel=r"$\cos\theta\ [-]$", ylabel="distribution [-]", xlim=(-1, 1), ylim=(0, 1)
    )

    x = np.array([float(s) for s in txt[0][1:-1]])
    pnum = 30
    #txtidx = np.linspace(1, len(txt)-1, num=pnum, dtype=int)
    txtidx = np.logspace(0, np.log10(len(txt)-1), num=pnum, dtype=int)
    for i, ti in enumerate(txtidx):
        timestep = float(txt[ti][0])
        y = np.array([float(s) for s in txt[ti][1:-1]])
        y /= simpson(y, x)

        # fitting
        popt, pcov = curve_fit(func, x, y, maxfev=10000, bounds=([0., 0, 2], [100, 1, 20]))
        y = func(x, *popt)
        strain = timestep * dt * PRM.fs2ns * erate
        ax.plot(
            x,
            y,
            label=r"$\varepsilon =$" + f"{strain}" + r"$[-]$",
            color=cm.jet((i - 1) / (len(txtidx) - 1)),
        )

    #ax.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=8)
    fig.savefig("sample_data/Cry.png", bbox_inches="tight")
