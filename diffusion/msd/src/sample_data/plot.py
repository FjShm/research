import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


if __name__ == "__main__":
    cg = pd.read_csv("msd.txt", sep="\s+")

    fig = plt.figure()
    ax = fig.add_subplot(
        111,
        xlabel=r"$time\ {\rm [ns]}$",
        ylabel=r"$MSD\ {\rm [\AA^2]}$",
        #xlim=(0.1, 500),
        #ylim=(1e-4, 1.2e0),
        #yticks=np.logspace(-4, 0, 5),
    )
    #tckr = ticker.LogLocator(base=10.0, subs="all", numticks=12)
    #ax.yaxis.set_minor_locator(tckr)
    ax.plot(cg["time"], cg["average"], label="CG")
    ax.legend()
    fig.savefig("compare.png", bbox_inches="tight")
