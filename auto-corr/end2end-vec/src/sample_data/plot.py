import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    cg = pd.read_csv("ete_auto-corr_fix.txt", header=None)
    fa = pd.read_csv("pure-cis-Ct-1per.dat", header=None, sep="\s+")
    tsf = float(sys.argv[1])
    fig = plt.figure(dpi=600)
    ax = fig.add_subplot(111,
            xlabel=r"unscaled $t\ {\rm [ns]}$",
            ylabel=r"$C(t)\ [-]$")
    ax.plot(cg[0]*tsf, cg[2], label="CG")
    ax.plot(fa[0], fa[1], label="FA")
    ax.legend()
    fig.savefig("test.png")
