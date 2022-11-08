import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    cg = pd.read_csv("rouse_mode-p_autocorr.txt", header=None, sep="\s+")
    fa = pd.read_csv("pure-cis-Ct-1per.dat", header=None, sep="\s+")
    fa[0] -= 8.70600033
    tsf = float(sys.argv[1])
    fig = plt.figure(dpi=600)
    ax = fig.add_subplot(
        111,
        xlabel=r"unscaled $t\ {\rm [ns]}$",
        ylabel=r"$C(t)\ [-]$",
        #xlim=(0, 400),
        xlim=(0.1, 500),
        ylim=(1e-5, 1e0),
        #yticks=np.arange(0, 1.1, 0.2),
        yticks=np.logspace(-5, 1, 7),
        xscale="log",
        yscale="log",
    )
    ax.plot(cg[0] * tsf, cg[1], label="CG")
    ax.plot(fa[0], fa[1], label="FA")
    ax.legend()
    fig.savefig("compare.png")
