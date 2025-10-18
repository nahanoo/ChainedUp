import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots

meta = pd.read_csv("../data/251018_succinate_glucose_plate_reader/metaod/metadata.csv")
data = pd.read_csv(
    "../data/251018_succinate_glucose_plate_reader/metaod/measurements.csv",
    index_col=0,
)
colors = {"Succinate": "#31688e", "Glucose": "#35b779", "Succinate+Glucose": "#fde725"}


def get_linegroups(species, cs, conc, meta):
    linegroups = meta[
        (meta["species"] == species)
        & (meta["carbon_source"] == cs)
        & (meta["cs_conc"] == conc)
    ]["linegroup"].unique()
    return linegroups


C_to_mM_succinate = {45: 11.25, 30: 7.5, 15: 3.75, 7.5: 1.875, 5: 1.25, 2.5: 0.625}
C_to_mM_glucose = {45: 7.5, 30: 5, 15: 2.5, 7.5: 1.25, 5: 0.833, 2.5: 0.417}
C_to_mM_glucose_succinate = {
    90: [11.25, 7.5],
    60: [7.5, 5],
    30: [3.75, 2.5],
    15: [1.875, 1.25],
    10: [1.25, 0.833],
    5: [0.625, 0.417],
}
species = ["At", "Oa", "At+Oa"]


def succinate():
    fig = make_subplots(rows=1, cols=3, subplot_titles=species, shared_xaxes=True)
    cs = "Succinate"
    for s in species:
        for c, conc in C_to_mM_succinate.items():
            linegroups = get_linegroups(s, cs, c, meta)
            for lg in linegroups:
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        line=dict(color="purple", width=2),
                        name=f"{conc} mM Succinate",
                        showlegend=False,
                        text=str(conc) + " mM",
                        textposition="middle center",
                    ),
                    row=1,
                    col=species.index(s) + 1,
                )
    fig.update_layout(
        title="Growth on Succinate",
        xaxis2_title="Time (hours)",
        yaxis_title="OD600",
    )
    fig = style_plot(fig)
    fig.write_image("plots/growth_curves/succinate.svg")


def glucose():
    fig = make_subplots(rows=1, cols=3, subplot_titles=species, shared_xaxes=True)
    cs = "Glucose"
    for s in species:
        for c, conc in C_to_mM_glucose.items():
            linegroups = get_linegroups(s, cs, c, meta)
            for lg in linegroups:
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        line=dict(color="purple", width=2),
                        name=f"{conc} mM Glucose",
                        showlegend=False,
                        text=str(conc) + " mM",
                        textposition="middle center",
                    ),
                    row=1,
                    col=species.index(s) + 1,
                )
    fig.update_layout(
        title="Growth on Glucose",
        xaxis2_title="Time (hours)",
        yaxis_title="OD600",
    )
    fig = style_plot(fig)
    fig.write_image("plots/growth_curves/glucose.svg")


def glucose_succinate():
    fig = make_subplots(rows=1, cols=3, subplot_titles=species, shared_xaxes=True)
    cs = "Succinate+Glucose"
    for s in species:
        for c, conc in C_to_mM_glucose_succinate.items():
            linegroups = get_linegroups(s, cs, c, meta)
            for lg in linegroups:
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        line=dict(color="purple", width=2),
                        name=f"{conc[0]} mM Succinate, {conc[1]} mM Glucose",
                        showlegend=False,
                        text=str(conc) + " mM",
                        textposition="middle center",
                    ),
                    row=1,
                    col=species.index(s) + 1,
                )
    fig.update_layout(
        title="Growth on Succinate+Glucose",
        xaxis2_title="Time (hours)",
        yaxis_title="OD600",
    )
    fig = style_plot(fig)
    fig.write_image("plots/growth_curves/glucose_succinate.svg")


fig = make_subplots(rows=1, cols=3, subplot_titles=species, shared_xaxes=True)
css = ["Glucose", "Succinate", "Succinate+Glucose"]
q = 30
concs = [q, q, 2 * q]
seen = set()
for i, cs in enumerate(css):
    for j, s in enumerate(species):
        linegroups = get_linegroups(s, cs, concs[i], meta)
        for lg in linegroups:
            x = data[lg + "_time"].to_numpy()
            y = data[lg + "_measurement"].to_numpy()
            show = cs not in seen
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line=dict(color=colors[cs], width=2),
                    name=cs,  # same name across traces of this group
                    legendgroup=cs,  # groups in legend
                    showlegend=show,  # only first one shows
                    # legendrank=10,         # optional: control ordering
                ),
                row=1,
                col=j + 1,
            )
            if show:
                seen.add(cs)
fig.update_layout(
    title="Growth on Succinate+Glucose",
    xaxis2_title="Time (hours)",
    yaxis_title="OD600",
)
fig = style_plot(fig)
fig.write_image("plots/growth_curves/glucose_succinate_layered.svg")
