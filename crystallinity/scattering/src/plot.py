import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from tqdm import tqdm


if __name__ == "__main__":
    files = glob.glob("results/data/joined/LOG.*.out")
    kx = np.loadtxt("results/data/kx.txt")
    ky = np.loadtxt("results/data/ky.txt").T
    interval = np.linspace(0, 10, 11)
    for f in tqdm(files):
        dat = np.loadtxt(f)
        fig = plt.figure()
        ax = fig.add_subplot(
            xlabel=r"$k_x\ {\rm [\AA^{-1}]}$",
            ylabel=r"$k_y\ {\rm [\AA^{-1}]}$",
            aspect="equal",
        )
        kx_, ky_ = np.meshgrid(kx, ky)
        cont = ax.contourf(kx_, ky_, dat, interval, cmap=cm.gray)
        fig.colorbar(cont)
        fname = os.path.basename(f).replace("LOG.", "").replace(".out", "")
        fig.savefig(f"results/{fname}.eps", bbox_inches="tight")
        fig.savefig(f"results/{fname}.png", bbox_inches="tight")
        plt.close(fig)
        del fig
