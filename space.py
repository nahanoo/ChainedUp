from scipy.integrate import odeint
import numpy as np
import plotly.graph_objects as go
import pandas as pd
from plotly.subplots import make_subplots
from style import *


def parse_params(cfus=False):
    if cfus:
        df = dict(pd.read_csv("parameters_cfus.csv"))
        params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
        p = params
        return p
    else:
        df = dict(pd.read_csv("parameters.csv"))
        params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
        p = params
        return p


def serial_c1c2(y, t, p):
    N1C1, R1C1, R2C1, N1C2, R1C2, R2C2 = y
    N1JC1 = p["v1_1"] * (R1C1 / (p["K1_1"] + R1C1))
    N1JC2 = p["v1_2"] * (R2C2 / (p["K1_2"] + R2C2))
    dN1C1 = N1C1 * (N1JC1 - p["D"])
    dR1C1 = -N1JC1 * N1C1 / p["q1_1"] + p["D"] * (p["M1"] - R1C1)
    dR2C1 = p["D"] * (p["M2"] - R2C1)
    dN1C2 = N1C2 * (N1JC2 - p["D"]) + p["D"] * N1C1
    dR2C2 = -N1JC2 * N1C2 / p["q1_2"] + p["D"] * (R2C1 - R2C2)
    dR1C2 = p["D"] * (R1C1 - R1C2)
    return dN1C1, dR1C1, dR2C1, dN1C2, dR1C2, dR2C2


xs = np.linspace(0, 200, 1000)
N1C1, R1C1, R2C1, N1C2, R1C2, R2C2 = odeint(
    serial_c1c2,
    [0.1, 10, 10, 0.1, 10, 10],
    xs,
    args=(parse_params(),),
).T

fig = go.Figure()
fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
for i, N in enumerate([N1C1, N1C2]):
    fig.add_trace(
        go.Scatter(x=xs, y=N, mode="lines", marker=dict(color="black")),
        row=1,
        col=1 + i,
    )
fig.update_layout(
    xaxis=dict(title="Time", ticks="inside"),
    yaxis=dict(title="OD<sub>600</sub>"),
    width=380,
    height=200,
    showlegend=False,
)
fig = style_plot(fig)
fig.write_image("plots/space_N1C1C2.svg")
