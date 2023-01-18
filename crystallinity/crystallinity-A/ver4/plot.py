import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pprint


class PRM:
    fs2ns = 1e-6
    chain_length = 49
    N = chain_length
    num_bond = chain_length - 1


if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], sep="\s+")
    erate = float(input("please type erate [ns^-1]: "))
    dt = float(input("type dt [fs]: "))

    mod_column_name = {}
    for h in df.columns:
        if h == "TimeStep":
            continue
        mod = int(h) % PRM.N
        if str(mod) in mod_column_name.keys():
            mod_column_name[str(mod)].append(h)
        else:
            mod_column_name[str(mod)] = [h]
        del mod

    # average columns of mod_column_name
    mean_col_names = []
    for key in mod_column_name.keys():
        dfc = df[key].copy()
        for col in mod_column_name[key][1:]:
            dfc += df[col]
        dfc /= len(mod_column_name[key])
        df[f"mean_{key}"] = dfc
        mean_col_names.append(f"mean_{key}")
        del dfc

    # distance from terminal
    distance_from_terminal = {}
    # {distance: ["col_name", "colname", ..]}
    for mcn in mean_col_names:
        idx = int(mcn.replace("mean_", ""))
        distance = min(idx, (PRM.N - 1) - idx)
        if distance in distance_from_terminal.keys():
            distance_from_terminal[distance].append(mcn)
        else:
            distance_from_terminal[distance] = [mcn]

    # prepare plot data
    timesteps = np.array(df["TimeStep"])
    strain = timesteps * dt * PRM.fs2ns * erate
    df = df[mean_col_names]
    x = [int(s.replace("mean_", "")) for s in df.columns]
    dist_from_trmnl_colname = []
    for key in distance_from_terminal.keys():
        dftc = f"dist_from_trmnl_{key}"
        df[dftc] = df[distance_from_terminal[key]].mean(axis=1)
        dist_from_trmnl_colname.append(dftc)
    df_dftc = df[dist_from_trmnl_colname]

    # plot
    fig = plt.figure(dpi=600)
    ax = fig.add_subplot(
        111,
        xlabel="strain [-]",
        ylabel="crystallinity [-]",
        xlim=(0, 14),
        xticks=np.arange(0, 15, 2),
        ylim=(0, 0.2),
        yticks=np.arange(0, 0.21, 0.05),
    )
    for i, col in enumerate(df_dftc.columns):
        y = list(df[col].rolling(10, center=True).mean())
        label = col.replace("dist_from_trmnl_", "")
        ax.plot(strain, y, label=label, color=cm.jet(i / len(df_dftc.columns)))
    ax.plot(-1, -1, color="white", label="# bonds from end")

    # legend関連設定
    hans, labs = ax.get_legend_handles_labels()
    ax.legend(
        bbox_to_anchor=(1, 1), loc="upper left", handles=hans[::-1], labels=labs[::-1]
    )
    Dir = os.path.dirname(sys.argv[1])
    fname = os.path.splitext(os.path.basename(sys.argv[1]))[0] + ".pdf"
    fig.savefig(os.path.join(Dir, fname), bbox_inches="tight")
