import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots
from fluorescence_od_mapping import map_od_to_gfp

meta = pd.read_csv("data/250918_glucose_lactose_plate/metadata.csv")
data = pd.read_csv("data/250918_glucose_lactose_plate/measurements.csv", index_col=0)


def get_linegroups(species, cs, signal, meta):
    linegroups = meta[
        (meta["species"] == species)
        & (meta["carbon_source"] == cs)
        & (meta["exp_ID"] == signal)
    ]["linegroup"].unique()
    return linegroups


colors = ["#8A2BE2", "#fc8d62", "#3D5BE2"]
height = 300
width = 300


def plot_at():
    species = "At"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_GFP",
    ]
    fig = make_subplots(
        cols=1,
        rows=2,
        vertical_spacing=0.05,
        shared_xaxes=True,
    )

    for i, signal in enumerate(signals):
        for j, cs in enumerate(css):
            lgs = get_linegroups(species, cs, signal, meta)
            for k, lg in enumerate(lgs):
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                if cs == "Glucose+Lactose":
                    y = 2 * y
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        name=cs,
                        legendgroup=cs,
                        marker=dict(color=colors[j]),
                        showlegend=False,
                    ),
                    row=i + 1,
                    col=1,
                )
    fig.update_layout(
        title_text="At growth and GFP expression",
        title_x=0.5,
        yaxis1_title="OD600",
        yaxis2_title="GFP",
        xaxis2_title="Time (h)",
        width=width,
        height=height,
    )
    fig = style_plot(fig)
    fig.write_image("plots/250918_glucose_lactose_at_od_gfp.svg")


def plot_oa():
    species = "Oa"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_mcherry",
    ]
    fig = make_subplots(
        cols=1,
        rows=2,
        vertical_spacing=0.05,
        shared_xaxes=True,
    )

    for i, signal in enumerate(signals):
        for j, cs in enumerate(css):
            lgs = get_linegroups(species, cs, signal, meta)
            for k, lg in enumerate(lgs):
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                if cs == "Glucose+Lactose":
                    y = 2 * y
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        name=cs,
                        legendgroup=cs,
                        marker=dict(color=colors[j]),
                        showlegend=False,
                    ),
                    row=i + 1,
                    col=1,
                )
    fig.update_layout(
        title_text="Oa growth and mCherry expression",
        title_x=0.5,
        yaxis1_title="OD600",
        yaxis2_title="mCherry",
        xaxis1_title="Time (h)",
        xaxis2_title="Time (h)",
        width=width,
        height=height,
    )
    fig = style_plot(fig)
    fig.write_image("plots/250918_glucose_lactose_oa_od_mcherry.svg")


def plot_at_oa():
    species = "At+Oa"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_mcherry",
        "glucose_lactose_250918_glucose_lactose_screening_GFP",
    ]
    fig = make_subplots(cols=1, rows=3, vertical_spacing=0.05)

    for i, signal in enumerate(signals):
        for j, cs in enumerate(css):
            lgs = get_linegroups(species, cs, signal, meta)
            for k, lg in enumerate(lgs):
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                if cs == "Glucose+Lactose":
                    y = 2 * y
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        name=cs,
                        legendgroup=cs,
                        marker=dict(color=colors[j]),
                        showlegend=False,
                    ),
                    row=i + 1,
                    col=1,
                )
    fig.update_layout(
        title_text="At+Oa",
        title_x=0.5,
        yaxis1_title="OD600",
        yaxis2_title="mCherry",
        yaxis3_title="GFP",
        xaxis3_title="Time (h)",
        width=width,
        height=3 / 2 * height,
    )
    fig = style_plot(fig)
    fig.write_image("plots/250918_glucose_lactose_at_oa_od_mcherry_gfp.svg")


plot_at()
plot_oa()
plot_at_oa()
species = "At+Ct+Ml+Oa"
css = ["Glucose", "Lactose", "Glucose+Lactose"]
signals = [
    "glucose_lactose_250918_glucose_lactose_screening_OD",
    "glucose_lactose_250918_glucose_lactose_screening_mcherry",
    "glucose_lactose_250918_glucose_lactose_screening_GFP",
]
fig = make_subplots(cols=1, rows=3, vertical_spacing=0.05)

for i, signal in enumerate(signals):
    for j, cs in enumerate(css):
        lgs = get_linegroups(species, cs, signal, meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()
            y = data[lg + "_measurement"].to_numpy()
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    name=cs,
                    legendgroup=cs,
                    marker=dict(color=colors[j]),
                    showlegend=False,
                ),
                row=i + 1,
                col=1,
            )
fig.update_layout(
    title_text="At+Ct+Oa",
    title_x=0.5,
    yaxis1_title="OD600",
    yaxis2_title="mCherry",
    yaxis3_title="GFP",
    xaxis3_title="Time (h)",
    width=width,
    height=3 / 2 * height,
)
fig = style_plot(fig)
fig.write_image("plots/250918_glucose_lactose_at_ct_oa_od_mcherry_gfp.svg")
