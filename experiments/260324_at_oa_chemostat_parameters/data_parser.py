from chibio_parser import *

reactor_mapping = {
    "M0": "Succinate",
    "M1": "Glucose",
    "M2": "Succinate+Glucose",
    "M3": "Succinate+Glucose Outflow",
}
species_mapping = {"at": "at_coculture", "oa": "oa_coculture"}

dfs = []
df, _ = cfu_parser("251208_at_chemostat_d_03")
dfs.append(df)
df, _ = cfu_parser("260119_oa_chemostat_d_03")
dfs.append(df)
df, _ = cfu_parser("260130_at_oa_chemostat_d_03")
df["species"] = df["species"].replace(species_mapping)
dfs.append(df)

df = pd.concat(dfs)
df["reactor"] = df["reactor"].replace(reactor_mapping)
df.to_csv("cfus.csv", index=False)
