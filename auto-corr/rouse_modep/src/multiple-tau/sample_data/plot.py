import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


if __name__ == "__main__":
    cg = pd.read_csv("rouse_mode-p_autocorr.txt", header=None, sep="\s+")
    fa = pd.read_csv("pure-cis-Ct-1per.dat", header=None, sep="\s+")
    oh = pd.read_csv("ohkuma_CGpc_rescaled.csv", header=None)

    fa[0] -= 8.70600033
    fig = plt.figure()
    ax = fig.add_subplot(
        111,
        xlabel=r"$st\ {\rm [ns]}$",
        ylabel=r"$C(t)\ [-]$",
        #xlim=(0, 400),
        xlim=(0.1, 500),
        ylim=(1e-4, 1.2e0),
        #yticks=np.arange(0, 1.1, 0.2),
        yticks=np.logspace(-4, 0, 5),
        xscale="log",
        yscale="log",
    )
    tckr = ticker.LogLocator(base=10.0, subs="all", numticks=12)
    ax.yaxis.set_minor_locator(tckr)
    tsf = float(sys.argv[1])
    ax.plot(cg[0] * tsf, cg[1], label=rf"AB-model ($s=${tsf})")
    ax.plot(fa[0], fa[1], label=r"FA ($s=$1)")
    ax.plot(oh[0], oh[1], label=r"CGpc by Ohkuma ($s=$8.3)")
    ax.legend()
    fig.savefig("compare.png", bbox_inches="tight")
