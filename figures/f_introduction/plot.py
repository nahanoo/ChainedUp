import pandas as pd
import plotly.graph_objects as go
from style import *
from model import Model
from experiment import Species
from os import path
import numpy as np


def plot_batch_growth_curves():
    meta = pd.read_csv(
        "../../data/251018_succinate_glucose_plate_reader/metaod/metadata.csv"
    )
    data = pd.read_csv(
        "../../data/251018_succinate_glucose_plate_reader/metaod/measurements.csv",
        index_col=0,
    )

    def get_linegroups(species, cs, conc, meta, signal):
        linegroups = meta[
            (meta["species"] == species)
            & (meta["carbon_source"] == cs)
            & (meta["cs_conc"] == conc)
            & (meta["exp_ID"] == signal)
        ]["linegroup"].unique()
        return linegroups

    carbon_sources = ["Succinate", "Glucose", "Succinate+Glucose"]
    cs_concs = [15, 15, 30]
    signals = [
        "2510_succinate_glucose_251018_succinate",
        "2510_succinate_glucose_251018_glucose",
        "2510_succinate_glucose_251018_succinate_glucose",
    ]
    fig = go.Figure()
    for i, carbon_source in enumerate(carbon_sources):
        linegroups = get_linegroups("At", carbon_source, cs_concs[i], meta, signals[i])
        for lg in linegroups:
            x = data[lg + "_time"].to_numpy()
            y = data[lg + "_measurement"].to_numpy()
            mask = x <= 27
            x = x[mask][::3]
            y = y[mask][::3]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    name=f"{carbon_source}",
                    line=dict(color=colors[carbon_source]),
                    showlegend=False,
                )
            )
    fig.update_layout(
        xaxis=dict(title="Time (h)"),
        yaxis=dict(title="OD<sub>600</sub>", range=[0, 0.8]),
        height=200,
        width=250,
        shapes=[
            dict(
                type="rect",
                yref="paper",
                xref="x",
                x0=0,
                x1=10,
                y0=0,
                y1=1,
                fillcolor=colors["Succinate"],
                opacity=0.3,
                layer="below",
                line_width=0,
            ),
            dict(
                type="rect",
                yref="paper",
                xref="x",
                x0=10,
                x1=max(x),
                y0=0,
                y1=1,
                fillcolor=colors["Glucose"],
                opacity=0.3,
                layer="below",
                line_width=0,
            ),
            dict(
                type="line",
                xref="x",
                yref="paper",
                x0=10,
                y0=0,
                x1=10,
                y1=1,
                line=dict(color="black", width=1, dash="dash"),
            ),
        ],
    )
    fig = style_plot(fig, line_thickness=2, right_margin=10, top_margin=10)
    fig.write_image("batch_growt_curves.svg")


conc = 7.5
p_f = path.join(
    "~", "ChainedUp", "modeling", "parameters", f"parameters_{conc}_mM_C.csv"
)
params = pd.read_csv(p_f, index_col=0)
figs = [go.Figure() for _ in range(3)]
a_params = [0, 1, 0.5]
fnames = ["allocation_0.svg", "allocation_1.svg", "allocation_0.5.svg"]
for i, a in enumerate(a_params):
    at, oa = Species("At", params.loc["At"]), Species("Oa", params.loc["Oa"])
    oa.N0 = 0
    at.a = a
    m = Model(at, oa, None, np.linspace(0, 100, 1000), 15, 0.15)
    m.integrate_model()
    figs[i].add_trace(
        go.Scatter(
            x=m.t,
            y=m.succinate,
            mode="lines",
            name="Succinate",
            line=dict(color=colors["Succinate"]),
        )
    )
    figs[i].add_trace(
        go.Scatter(
            x=m.t,
            y=m.glucose,
            mode="lines",
            name="Glucose",
            line=dict(color=colors["Glucose"]),
        )
    )
    figs[i].update_layout(
        xaxis=dict(title="Time (h)"),
        yaxis=dict(title=""),
        height=110,
        width=110,
        showlegend=False,
    )
    figs[i] = style_plot(
        figs[i], line_thickness=2, left_margin=0, right_margin=10, top_margin=10
    )
    figs[i].write_image(fnames[i])
fig = go.Figure()
at, oa = Species("At", params.loc["At"]), Species("Oa", params.loc["Oa"])
oa.N0 = 0
at.a = 0.5
m = Model(at, oa, None, np.linspace(0, 100, 1000), 15, 0.15)
m.integrate_model()
fig.add_trace(
    go.Scatter(
        x=m.t,
        y=m.at.y,
        mode="lines",
        name="Succinate",
        line=dict(color=colors["At"]),
    )
)
fig.update_layout(
    xaxis=dict(title="Time (h)"),
    yaxis=dict(title="OD<sub>600</sub>"),
    height=110,
    width=110,
    showlegend=False,
)
fig = style_plot(fig, line_thickness=2, top_margin=10)
fig.write_image("allocation_0.5_at.svg")
