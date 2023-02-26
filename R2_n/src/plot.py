import os
import sys
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("cannot read input file.")
        exit()
    dat = np.loadtxt(sys.argv[1])
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(
        111,
        xlabel=r"$n$",
        ylabel=r"$\langle R(n)^2\rangle/n$",
        xlim=(0, 50),
        ylim=(0.06, 0.2),
        xticks=np.arange(0, 51, 5),
        yticks=np.arange(0.06, 0.21, 0.02),
    )
    ax.plot(dat[:,0], dat[:,1]*0.01, "r-")
    fig.savefig(f"{os.path.splitext(sys.argv[1])[0]}.eps", bbox_inches="tight")
