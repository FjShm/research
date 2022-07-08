import pandas as pd

names = ("AA", "AB", "BB")
for name in names:
    df = pd.read_table(f"rdf.{name}", sep="\s+")
    df["r"] = df["r"] * 10
    df.to_csv(f"rdf.{name}", index=None, sep="\t")
