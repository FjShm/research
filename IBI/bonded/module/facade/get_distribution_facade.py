import pandas as pd


class GetDistributionFacade:
    def __call__(self, dist_adapter):
        path = dist_adapter.path
        ratio = dist_adapter.ratio
        df = pd.read_table(path, header=None, skiprows=1, sep="\s+")
        thre = df[1].max() * ratio
        df = df[df[1] >= thre]
        x = list(df[0])
        P = list(df[1])
        return x, P
