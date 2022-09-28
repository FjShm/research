import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("cannot read input file.")
        eixt()
    df = pd.read_csv(sys.argv[1], sep="\s+", header=None)
    fig = plt.figure(dpi=500)
    ax = fig.add_subplot(
        111,
        xlabel=r"$n$",
        ylabel=r"$\langle R(n)^2\rangle/n$",
        xlim=(0, 50),
        ylim=(0.06, 0.2),
        xticks=np.arange(0, 51, 5),
        yticks=np.arange(0.06, 0.21, 0.02),
    )
    ax.plot(list(df[0]), list(df[1] / 100), "r-")
    fig.savefig(os.path.splitext(sys.argv[1])[0], bbox_inches="tight")
