import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    cg = pd.read_csv("ete_auto-corr.txt", header=None, sep="\s+")
    fa = pd.read_csv("pure-cis-Ct-1per.dat", header=None, sep="\s+")
    fig = plt.figure()
    ax = fig.add_subplot(111,
            xlabel=r"$st\ {\rm [ns]}$",
            ylabel=r"$C(t)=\frac{\langle R(t) R(0)\rangle}{\langle R(0)^2\rangle}\ [-]$",
                         )
    tsf = float(sys.argv[1])
    ax.plot(cg[0]*tsf, cg[2], label=rf"AB-model ($s=${tsf})")
    ax.plot(fa[0], fa[1], label=r"FA ($s=$1)")
    ax.legend()
    fig.savefig("compare.png", bbox_inches="tight")
    fig.savefig("compare.eps", bbox_inches="tight")
