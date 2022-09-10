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
    df.sort_values(by=len(df)-1, axis=1, ascending=False, inplace=True)
    erate = float(input("please type erate [ns^-1]: "))
    dt = float(input("type dt [fs]: "))

    fig = plt.figure(dpi=600)
    ax = fig.add_subplot(111, xlabel=r"strain $\varepsilon$", ylabel="crystallinity [-]",
            xscale="log")
    for i, col in enumerate(df.columns):
        if i == 0:
            continue
        k = int((PRM.num_bond - (len(df.columns) - 1)) * 0.5)
        idx = k + col - 1
        distance_from_terminal = min(idx, PRM.num_bond - 1 - idx)
        ax.plot(
            list(df[0] * dt * PRM.fs2ns * erate),
            list(df[col]),
            label=f"position-{distance_from_terminal}",
            color=cm.jet((i-1)/(len(df.columns)-1)),
        )
    ax.legend()
    fig.savefig("sample_data/Cry.png", bbox_inches="tight")
