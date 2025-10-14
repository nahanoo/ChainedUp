import pandas as pd
from glob import glob
import plotly.graph_objects as go
import os
import re

time_points = {}
fs = [os.path.split(f)[-1] for f in glob("OD600_wholePlate/P2*.csv")]
for f in fs:
    t = f.split("-")[0]
    time_points[f] = int(t[t.find("T") + 1 :])

time_points = dict(sorted(time_points.items(), key=lambda x: x[1]))
A7 = []
A8 = []
A9 = []
for key in time_points.keys():
    df = pd.read_csv("OD600_wholePlate/" + key, index_col=0)
    A7.append(df.at["A", " 7"])
    A8.append(df.at["A", " 8"])
    A9.append(df.at["A", " 9"])


fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=list(time_points.values()),
        y=A7,
        mode="lines+markers",
        name="A7",
        line=dict(color="#1f77b4"),
    )
)
fig.add_trace(
    go.Scatter(
        x=list(time_points.values()),
        y=A8,
        mode="lines+markers",
        name="A8",
        line=dict(color="#ff7f0e"),
    )
)
fig.add_trace(
    go.Scatter(
        x=list(time_points.values()),
        y=A9,
        mode="lines+markers",
        name="A9",
        line=dict(color="#2ca02c"),
    )
)
fig.show()
