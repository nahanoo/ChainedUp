from chibio_parser import *

mapping = {
    "M0": "Succinate",
    "M1": "Glucose",
    "M2": "Succinate+Glucose",
    "M3": "Succinate+Glucose Outflow",
}

dfs = []
df, _ = cfu_parser("251208_at_chemostat_d_03")
dfs.append(df)
df, _ = cfu_parser("260119_oa_chemostat_d_03")
dfs.append(df)

df = pd.concat(dfs)
df["reactor"] = df["reactor"].replace(mapping)
df.to_csv("cfus.csv", index=False)
