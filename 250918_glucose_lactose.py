import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots
from fluorescence_od_mapping import map_od_to_gfp

meta = pd.read_csv("data/250918_glucose_lactose_plate/metadata.csv")
data = pd.read_csv("data/250918_glucose_lactose_plate/measurements.csv", index_col=0)

df_mcherry = pd.read_csv("data/mcherry_to_od.csv")


def mcherry_to_od(cs, values):
    map_od_to_gfp(cs)
    df_mcherry = pd.read_csv("data/mcherry_to_od.csv")
    mcherry_lookup = df_mcherry["mcherry"].values
    od_from_mcherry = df_mcherry["od"].values
    return np.interp(values, mcherry_lookup, od_from_mcherry)


def get_linegroups(species, cs, signal, meta):
    linegroups = meta[
        (meta["species"] == species)
        & (meta["carbon_source"] == cs)
        & (meta["exp_ID"] == signal)
    ]["linegroup"].unique()
    return linegroups


def reconstruct_od_from_gfp(x, gfp_values, od_init):
    x_keep = []
    y_keep = []
    for xi, yi in zip(x, gfp_values):
        if yi > 0:
            x_keep.append(xi)
            y_keep.append(yi)
    N_reconstructed = [od_init]
    r_gfp = np.diff(np.log(y_keep)) / np.diff(x_keep) / 1.4
    for i, m in enumerate(r_gfp[:-1]):
        dt = x_keep[i + 1] - x_keep[i]
        N_next = N_reconstructed[-1] * np.exp(m * dt)
        N_reconstructed.append(N_next)
    N_reconstructed = np.array(N_reconstructed)
    return N_reconstructed


def plot_at():
    species = "At"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_GFP",
    ]
    colors_array = ["black", "black"]

    fig = make_subplots(
        rows=2,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        horizontal_spacing=0.05,
        vertical_spacing=0.1,
        subplot_titles=css,
    )
    for j, cs in enumerate(css):
        lgs = get_linegroups(species, cs, signals[0], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            od_init = y[0]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="blue"),
                    name="OD600",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[1], meta)
        y_gfps = []
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="GFP",
                    showlegend=False,
                ),
                row=2,
                col=j + 1,
            )
            y_gfps.append(y)
        y_gfp = np.mean(y_gfps, axis=0)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=reconstruct_od_from_gfp(x, y_gfp, od_init),
                mode="lines",
                marker=dict(color=colors["at"]),
                name="At (OD<br>from GFP)",
                showlegend=(True if j + k == 0 else False),
            ),
            row=1,
            col=j + 1,
        )
    fig.update_layout(
        xaxis5=dict(title="Time (hours)"),
        yaxis1=dict(title="OD600"),
        yaxis4=dict(title="GFP"),
        title="Glucose/Lactose screening of At",
    )
    fig = style_plot(fig, font_size=10, marker_size=4)
    fig.write_image("plots/250918_glucose_lactose_at.svg")


def plot_oa():
    species = "Oa"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_mcherry",
    ]

    fig = make_subplots(
        rows=2,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        horizontal_spacing=0.05,
        vertical_spacing=0.1,
        subplot_titles=css,
    )
    for j, cs in enumerate(css):
        lgs = get_linegroups(species, cs, signals[0], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="blue"),
                    name="OD600",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[1], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="mcherry",
                    showlegend=False,
                ),
                row=2,
                col=j + 1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mcherry_to_od(cs, y),
                    mode="lines",
                    marker=dict(color=colors["oa"]),
                    name="Oa (OD<br>from mcherry)",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
    fig.update_layout(
        xaxis5=dict(title="Time (hours)"),
        yaxis1=dict(title="OD600"),
        yaxis4=dict(title="GFP"),
        title="Glucose/Lactose screening of Oa",
    )
    fig = style_plot(fig, font_size=10, marker_size=4)
    fig.write_image("plots/250918_glucose_lactose_oa.svg")


def plot_at_oa():
    species = "At+Oa"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_mcherry",
        "glucose_lactose_250918_glucose_lactose_screening_GFP",
    ]

    fig = make_subplots(
        rows=3,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        horizontal_spacing=0.05,
        vertical_spacing=0.1,
        subplot_titles=css,
    )
    for j, cs in enumerate(css):
        lgs = get_linegroups(species, cs, signals[0], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            od_init = y[0]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="blue"),
                    name="OD600",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[1], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="mcherry",
                    showlegend=False,
                ),
                row=2,
                col=j + 1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mcherry_to_od(cs, y),
                    mode="lines",
                    marker=dict(color=colors["oa"]),
                    name="Oa (OD<br>from mcherry)",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[2], meta)
        y_gfps = []
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="GFP",
                    showlegend=False,
                ),
                row=2,
                col=j + 1,
            )
            y_gfps.append(y)
        y_gfp = np.mean(y_gfps, axis=0)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=reconstruct_od_from_gfp(x, y_gfp, od_init),
                mode="lines",
                marker=dict(color=colors["at"]),
                name="At (OD<br>from GFP)",
                showlegend=(True if j + k == 0 else False),
            ),
            row=1,
            col=j + 1,
        )
    fig.update_layout(
        xaxis5=dict(title="Time (hours)"),
        yaxis1=dict(title="OD600"),
        yaxis4=dict(title="mcherry"),
        yaxis7=dict(title="GFP"),
        title="Glucose/Lactose screening of At+Oa",
    )
    fig = style_plot(fig, font_size=10, marker_size=4)
    fig.write_image("plots/250918_glucose_lactose_at_oa.svg")


plot_at_oa()


def plot_at_ml_oa():
    species = "At+Ml+Oa"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_mcherry",
        "glucose_lactose_250918_glucose_lactose_screening_GFP",
    ]

    fig = make_subplots(
        rows=3,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        horizontal_spacing=0.05,
        vertical_spacing=0.1,
        subplot_titles=css,
    )
    for j, cs in enumerate(css):
        lgs = get_linegroups(species, cs, signals[0], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="blue"),
                    name="OD600",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[1], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="mcherry",
                    showlegend=False,
                ),
                row=2,
                col=j + 1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mcherry_to_od(cs, y),
                    mode="lines",
                    marker=dict(color=colors["oa"]),
                    name="Oa (OD<br>from mcherry)",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[2], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="GFP",
                    showlegend=False,
                ),
                row=3,
                col=j + 1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=gfp_to_od(cs, y),
                    mode="lines",
                    marker=dict(color=colors["at"]),
                    name="At (OD<br>from GFP)",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
    fig.update_layout(
        xaxis5=dict(title="Time (hours)"),
        yaxis1=dict(title="OD600"),
        yaxis4=dict(title="mcherry"),
        yaxis7=dict(title="GFP"),
        title="Glucose/Lactose screening of At+Ml+Oa",
    )

    fig = style_plot(fig, font_size=10, marker_size=4)
    fig.write_image("plots/250918_glucose_lactose_at_ml_oa.svg")


def plot_at_ct_ml_oa():
    species = "At+Ct+Ml+Oa"
    css = ["Glucose", "Lactose", "Glucose+Lactose"]
    signals = [
        "glucose_lactose_250918_glucose_lactose_screening_OD",
        "glucose_lactose_250918_glucose_lactose_screening_mcherry",
        "glucose_lactose_250918_glucose_lactose_screening_GFP",
    ]

    fig = make_subplots(
        rows=3,
        cols=3,
        shared_yaxes=True,
        shared_xaxes=True,
        horizontal_spacing=0.05,
        vertical_spacing=0.1,
        subplot_titles=css,
    )
    for j, cs in enumerate(css):
        lgs = get_linegroups(species, cs, signals[0], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="blue"),
                    name="OD600",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[1], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="mcherry",
                    showlegend=False,
                ),
                row=2,
                col=j + 1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mcherry_to_od(cs, y),
                    mode="lines",
                    marker=dict(color=colors["oa"]),
                    name="Oa (OD<br>from mcherry)",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
        lgs = get_linegroups(species, cs, signals[2], meta)
        for k, lg in enumerate(lgs):
            x = data[lg + "_time"].to_numpy()[1:]
            y = data[lg + "_measurement"].to_numpy()[1:]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    marker=dict(color="black"),
                    name="GFP",
                    showlegend=False,
                ),
                row=3,
                col=j + 1,
            )
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=gfp_to_od(cs, y),
                    mode="lines",
                    marker=dict(color=colors["at"]),
                    name="At (OD<br>from GFP)",
                    showlegend=(True if j + k == 0 else False),
                ),
                row=1,
                col=j + 1,
            )
    fig.update_layout(
        xaxis5=dict(title="Time (hours)"),
        yaxis1=dict(title="OD600"),
        yaxis4=dict(title="mcherry"),
        yaxis7=dict(title="GFP"),
        title="Glucose/Lactose screening of At+Ct+Ml+Oa",
    )
    fig = style_plot(fig, font_size=10, marker_size=4)
    fig.write_image("plots/250918_glucose_lactose_at_ct_ml_oa.svg")
