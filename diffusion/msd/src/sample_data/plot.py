import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

if __name__ == "__main__":
    cg = pd.read_csv("msd.txt", sep="\s+")

    fig = plt.figure()
    ax = fig.add_subplot(
        111,
        xlabel=r"$time\ {\rm [ns]}$",
        ylabel=r"$MSD\ {\rm [\AA^2]}$",
        #xlim=(0.1, 500),
        #ylim=(1e-4, 1.2e0),
    )
    cols = [col for col in cg.columns if col not in ("time", "average")]
    for i,col in enumerate(cols):
        ax.plot(cg["time"], cg[col], color=cm.jet(i/(len(cols)-1)), linewidth=1)
    ax.plot(cg["time"], cg["average"], color="k")
    #ax.legend()
    fig.savefig("compare.png", bbox_inches="tight")
