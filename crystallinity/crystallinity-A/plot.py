import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class PRM:
    fs2ns = 1e-6
    chain_length = 49
    num_bond = chain_length - 1


if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], header=None, sep="\s+")
    df.sort_values(by=len(df) - 1, axis=1, ascending=False, inplace=True)
    erate = float(input("please type erate [ns^-1]: "))
    dt = float(input("type dt [fs]: "))

    fig = plt.figure(dpi=600)
    ax = fig.add_subplot(
        111, xlabel=r"strain $\varepsilon$", ylabel="crystallinity [-]", xscale="log"
    )
    ax.plot([1], [0.8], color="white", marker=".", label="# bonds from terminal")
    for i, col in enumerate(df.columns):
        if i == 0:
            continue
        k = int((PRM.num_bond - (len(df.columns) - 1)) * 0.5)
        idx = k + col - 1
        distance_from_terminal = min(idx, PRM.num_bond - 1 - idx)
        ax.plot(
            list(df[0] * dt * PRM.fs2ns * erate),
            list(df[col]),
            label=f"{distance_from_terminal:^10d}",
            color=cm.jet((i - 1) / (len(df.columns) - 1)),
        )
    ax.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=12)
    Dir = os.path.dirname(sys.argv[1])
    fig.savefig(os.path.join(Dir, "Cry.png"), bbox_inches="tight")
