from chibio_parser import cfu_parser
import pandas as pd

dfs = []
df, _ = cfu_parser("251208_at_chemostat_d_03")
df["comment"] = "At"
dfs.append(df)
df, _ = cfu_parser("260119_oa_chemostat_d_03")
df["comment"] = "Oa"
dfs.append(df)
df, _ = cfu_parser("260130_at_oa_chemostat_d_03")
df["comment"] = "AtOA"
dfs.append(df)
pd.concat(dfs).to_csv("data.csv", index=False)
