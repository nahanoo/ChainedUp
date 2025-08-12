import numpy as np
from scipy.integrate import odeint
from style import *
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


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


def N1C1C2(y, t, p):
    N1C1, N1C2, R1C1, R1C2 = y
    N1JC1 = p["v1_1"] * R1C1 / (R1C1 + p["K1_1"])
    N1JC2 = p["v1_1"] * R1C2 / (R1C2 + p["K1_1"])
    dN1C1 = N1C1 * (N1JC1 - p["D"])
    dN1C2 = N1C2 * (N1JC2 - p["D"]) + p["D"] * N1C1
    dR1C1 = p["D"] * (p["M1"] - R1C1) - N1JC1 * N1C1 / p["q1_1"]
    dR1C2 = p["D"] * (R1C1 - R1C2) - N1JC2 * N1C2 / p["q1_1"]
    return [dN1C1, dN1C2, dR1C1, dR1C2]


def N1C1C2_strategy(y, t, p):
    N1C1, N1C2, R1C1, R1C2, M1C1, M1C2 = y
    N1JC1 = p["v1_1"] * R1C1 / (R1C1 + p["K1_1"])
    N1JC2 = p["v1_1"] * R1C2 / (R1C2 + p["K1_1"])
    dN1C1 = N1C1 * (N1JC1 - p["D"])
    dN1C2 = (
        N1C2 * (N1JC2 - p["D"]) + p["D"] * N1C1 + p["v1_1"] * M1C2 / (M1C2 + p["K1_2"])
    )
    dR1C1 = p["D"] * (p["M1"] - R1C1) - N1JC1 * N1C1 / p["q1_1"]
    dR1C2 = p["D"] * (R1C1 - R1C2) - N1JC2 * N1C2 / p["q1_1"]
    dM1C1 = N1JC1 * N1C1 / p["a1_2"] - p["D"] * M1C1
    dM1C2 = (
        p["D"] * (M1C1 - M1C2)
        - p["v1_1"] * M1C2 / (M1C2 + p["K1_2"]) * N1C2 / p["q1_2"]
    )
    return [dN1C1, dN1C2, dR1C1, dR1C2, dM1C1, dM1C2]


def plot_N1C1C2():
    p = parse_params()
    xs = np.linspace(0, 1000, 10000)
    y0 = [1, 0, 10, 10]
    ys = odeint(N1C1C2, y0, xs, args=(p,))
    N1C1, N1C2, R1C1, R1C2 = ys.T
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
    fig.write_image("plots/N1C1C2.svg")
    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    for i, R in enumerate([R1C1, R1C2]):
        fig.add_trace(
            go.Scatter(x=xs, y=R, mode="lines", marker=dict(color="black")),
            row=1,
            col=1 + i,
        )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="Concentration [mM]", type="log"),
        width=380,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/R1C1C2.svg")
    print("N1C1:", N1C1[-1], "N1C2:", N1C2[-1])
    print("R1C1:", R1C1[-1], "R1C2:", R1C2[-1])


def N1N2C1C2(y, t, p):
    N1C1, N2C1, N1C2, N2C2, R1C1, R1C2 = y
    N1JC1 = p["v1_1"] * R1C1 / (R1C1 + p["K1_1"])
    N1JC2 = p["v1_1"] * R1C2 / (R1C2 + p["K1_1"])
    N2JC1 = p["v2_1"] * R1C1 / (R1C1 + p["K2_1"])
    N2JC2 = p["v2_1"] * R1C2 / (R1C2 + p["K2_1"])
    dN1C1 = N1C1 * (N1JC1 - p["D"])
    dN1C2 = N1C2 * (N1JC2 - p["D"]) + p["D"] * N1C1
    dN2C1 = N2C1 * (N2JC1 - p["D"])
    dN2C2 = N2C2 * (N2JC2 - p["D"]) + p["D"] * N2C1
    dR1C1 = (
        p["D"] * (p["M1"] - R1C1) - N1JC1 * N1C1 / p["q1_1"] - N2JC1 * N2C1 / p["q2_1"]
    )
    dR1C2 = p["D"] * (R1C1 - R1C2) - N1JC2 * N1C2 / p["q1_1"] - N2JC2 * N2C2 / p["q2_1"]
    return [dN1C1, dN2C1, dN1C2, dN2C2, dR1C1, dR1C2]


def plot_N1N2C1C2():
    p = parse_params()
    xs = np.linspace(0, 1000, 10000)
    y0 = [1, 1, 0, 0, 10, 10]
    ys = odeint(N1N2C1C2, y0, xs, args=(p,))
    N1C1, N2C1, N1C2, N2C2, R1C1, R1C2 = ys.T
    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    fig.add_trace(
        go.Scatter(x=xs, y=N1C1, mode="lines", marker=dict(color=colors["ct"])),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(x=xs, y=N2C1, mode="lines", marker=dict(color=colors["oa"])),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(x=xs, y=N1C2, mode="lines", marker=dict(color=colors["ct"])),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(x=xs, y=N2C2, mode="lines", marker=dict(color=colors["oa"])),
        row=1,
        col=2,
    )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="OD<sub>600</sub>"),
        width=380,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/N1N2C1C2.svg")


def plot_N1C1C2_cfus():
    p = parse_params(cfus=True)
    xs = np.linspace(0, 1000, 10000)
    y0 = [100000, 0, 10, 10]
    ys = odeint(N1C1C2, y0, xs, args=(p,))
    N1C1, N1C2, R1C1, R1C2 = ys.T
    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    for i, N in enumerate([N1C1, N1C2]):
        fig.add_trace(
            go.Scatter(x=xs, y=N, mode="lines", marker=dict(color="black")),
            row=1,
            col=1 + i,
        )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="CFUs/mL"),
        width=380,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/N1C1C2_cfus.svg")
    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    for i, R in enumerate([R1C1, R1C2]):
        fig.add_trace(
            go.Scatter(x=xs, y=R, mode="lines", marker=dict(color="black")),
            row=1,
            col=1 + i,
        )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="Concentration"),
        width=380,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/R1C1C2_cfus.svg")
    print("N1C1:", N1C1[-1], "N1C2:", N1C2[-1])
    print("R1C1:", R1C1[-1], "R1C2:", R1C2[-1])


def N1O2C1C2(y, t, p):
    N1C1, N1C2, R1C1, R1C2, OC1, OC2 = y
    N1R1JC1 = p["v1_1"] * R1C1 / (R1C1 + p["K1_1"])
    N1R1JC2 = p["v1_1"] * R1C2 / (R1C2 + p["K1_1"])
    N1OJC1 = OC1 / (OC1 + p["K1_O"])
    N1OJC2 = OC2 / (OC2 + p["K1_O"])
    N1JC1 = N1R1JC1 * N1OJC1
    N1JC2 = N1R1JC2 * N1OJC2
    dN1C1 = N1C1 * (N1JC1 - p["D"])
    dN1C2 = N1C2 * (N1JC2 - p["D"]) + p["D"] * N1C1
    dR1C1 = p["D"] * (p["M1"] - R1C1) - N1JC1 * N1C1 / p["q1_1"]
    dR1C2 = p["D"] * (R1C1 - R1C2) - N1JC2 * N1C2 / p["q1_1"]
    dOC1 = p["kla_1"] * (p["osat_1"] - OC1) - N1JC1 * OC1 / p["q1_O"] - p["D"] * OC1
    dOC2 = p["kla_1"] * (p["osat_1"] - OC2) - N1JC2 * OC2 / p["q1_O"] - p["D"] * OC2
    return [dN1C1, dN1C2, dR1C1, dR1C2, dOC1, dOC2]


def plot_N1O2C1C2():
    p = parse_params()
    xs = np.linspace(0, 1000, 10000)
    y0 = [0.1, 0, 10, 10, 8, 8]
    ys = odeint(N1O2C1C2, y0, xs, args=(p,))
    N1C1, N1C2, R1C1, R1C2, OC1, OC2 = ys.T
    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    fig.add_trace(
        go.Scatter(x=xs, y=N1C1, mode="lines", marker=dict(color="black")),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(x=xs, y=N1C2, mode="lines", marker=dict(color="black")),
        row=1,
        col=2,
    )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="OD<sub>600</sub>"),
        width=380,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/N1O2C1C2.svg")
    print("R1C1:", R1C1[-1], "R1C2:", R1C2[-1], " OC1:", OC1[-1], "OC2:", OC2[-1])

    fig = go.Figure()

    JR = [p["v1_1"] * R / (R + p["K1_1"]) for R in R1C1]
    JO = [p["v1_1"] * O / (O + p["K1_O"]) for O in OC1]
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=JR,
            mode="lines",
            marker=dict(color=colors["ct"]),
            name="Resource<br>growth<br>rate",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=JO,
            mode="lines",
            marker=dict(color=colors["oa"]),
            name="Oxygen<br>growth<br>rate",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="Growth rate"),
        width=380 / 2,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/JRJO.svg")
    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    for i, R in enumerate([R1C1, R1C2]):
        fig.add_trace(
            go.Scatter(x=xs, y=R, mode="lines", marker=dict(color="black")),
            row=1,
            col=1 + i,
        )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="Concentration [mM]"),
        width=380,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/R1O2C1C2.svg")


def N1C1C2_cross_feeding():
    p = parse_params(cfus=True)
    p["D"] = 0.1
    xs = np.linspace(0, 200, 10000)
    y0 = [200000, 0, 10, 10, 0, 0]
    ys = odeint(N1C1C2_strategy, y0, xs, args=(p,))
    N1C1, N1C2, R1C1, R1C2, M1C1, M1C2 = ys.T
    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    for i, N in enumerate([N1C1, N1C2]):
        fig.add_trace(
            go.Scatter(x=xs, y=N, mode="lines", marker=dict(color="black")),
            row=1,
            col=1 + i,
        )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="CFUs/mL", type="log"),
        width=380,
        height=200,
        showlegend=False,
    )
    # fig.show()
    fig = style_plot(fig)
    fig.write_image("plots/N1C1C2_cross_feeding_D" + str(p["D"]) + ".svg")

    fig = make_subplots(cols=2, rows=1, subplot_titles=["C1", "C2"], shared_yaxes=True)
    for i, R in enumerate([M1C1, M1C2]):
        fig.add_trace(
            go.Scatter(x=xs, y=R, mode="lines", marker=dict(color="black")),
            row=1,
            col=1 + i,
        )
    fig.update_layout(
        xaxis=dict(title="Time", ticks="inside"),
        yaxis=dict(title="Concentration [mM]"),
        width=380,
        height=200,
        showlegend=False,
    )
    fig = style_plot(fig)
    fig.write_image("plots/M1C1C2_cross_feeding_D" + str(p["D"]) + ".svg")

    print("N1C1:", N1C1[-1], "N1C2:", N1C2[-1])
    print("R1C1:", R1C1[-1], "R1C2:", R1C2[-1])
    print("M1C1:", M1C1[-1], "M1C2:", M1C2[-1])
