import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm


if __name__ == "__main__":
    txt = np.loadtxt(sys.argv[1], dtype=float, delimiter=" ")
    if len(txt.shape) == 1:
        txt = [txt]
    for row in tqdm(txt):
        dim = int(np.sqrt(len(row) - 1))
        timestep, row2d = row[0], row[1:].reshape(dim, dim)
        fig = plt.figure(dpi=600)
        sns.heatmap(row2d, vmax=10, vmin=0)
        fig.savefig(f"sample_data/s_{timestep}.png", bbox_inches="tight")
        plt.close(fig)
