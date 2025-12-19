from model_chain import Model
from experiment import Species
import numpy as np
from style import *
import plotly.graph_objects as go
import pandas as pd
from os import path
from multiprocessing import Pool
from plotly.subplots import make_subplots

font_size = 14


def simulate_endpoints(args):
    conc, params, xs, species_name, a = args

    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])

    if species_name == "At":
        oa.N0 = 0.0
        at.a = a
    elif species_name == "Oa":
        at.N0 = 0.0
        oa.a = a
    elif species_name == "At+Oa":
        at.a, oa.a = a

    model = Model(at, oa, None, xs, conc, 0.3)
    model.integrate_model()

    return (
        model.c1.succinate[-1],
        model.c1.glucose[-1],
        model.c1.at_y[-1],
        model.c1.oa_y[-1],
        model.c2.succinate[-1],
        model.c2.glucose[-1],
        model.c2.at_y[-1],
        model.c2.oa_y[-1],
    )


def a_mono_culture():
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    xs = np.linspace(0, 100, 1000)
    aas = np.linspace(0.0, 1.0, 50)
    with Pool() as pool:
        args_at = [(conc, params, xs, "At", a) for a in aas]
        results_at = pool.map(simulate_endpoints, args_at)
        suc_at_c1, gluc_at_c1, at_c1, oa_c1, suc_at_c2, gluc_at_c2, at_c2, oa_c2 = (
            np.asarray(results_at).T
        )
        args_oa = [(conc, params, xs, "Oa", a) for a in aas]
        results_oa = pool.map(simulate_endpoints, args_oa)
        suc_oa_c1, gluc_oa_c1, _, oa_c1, suc_oa_c2, gluc_oa_c2, _, oa_c2 = np.asarray(
            results_oa
        ).T

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=suc_at_c1,
            y=gluc_at_c1,
            mode="markers",
            marker=dict(
                symbol="triangle-up",
                color=aas,
                colorscale="cividis",
                colorbar=dict(
                    title="a",
                    len=0.55,  # height as fraction of plotting area
                    lenmode="fraction",
                    y=0.5,
                    yanchor="middle",
                    thickness=14,  # width in px (optional)
                ),
            ),
            name="At",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=suc_oa_c1,
            y=gluc_oa_c1,
            mode="markers",
            marker=dict(
                color=aas,
                colorscale="cividis",
                showscale=False,
            ),
            line=dict(width=2),
            name="Oa",
        )
    )
    fig.update_layout(
        xaxis=dict(
            title="Succinate [mM]",
            # type="log"
        ),
        yaxis=dict(
            title="Glucose [mM]",
            # type="log"
        ),
        title="Mono-culture",
    )
    fig = style_plot(fig, line_thickness=3, marker_size=10, font_size=font_size)
    fig.write_image("plots/contours/mono_culture_resource_allocation.svg")


def resource_allocation_heatmap():
    colors_heatmap = [
        colors["Glucose"],
        colors["Succinate"],
    ]
    aas = np.linspace(0.0, 1.0, 20)
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    xs = np.linspace(0, 100, 1000)
    with Pool() as pool:
        args = [
            (conc, params, xs, "At+Oa", (a_at, a_oa)) for a_at in aas for a_oa in aas
        ]
        results = pool.map(simulate_endpoints, args)
        suc_c1, gluc_c1, at_c1, oa_c1, suc_c2, gluc_c2, at_c2, oa_c2 = np.asarray(
            results
        ).T

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("First chemostat", "Downstream chemostat"),
        horizontal_spacing=0.01,
        shared_xaxes=True,
        shared_yaxes=True,
    )
    n = len(aas)
    at_c1_grid, oa_c1_grid = at_c1.reshape(n, n).T, oa_c1.reshape(n, n).T
    z_c1 = np.log10(at_c1_grid / oa_c1_grid)
    at_c2_grid, oa_c2_grid = at_c2.reshape(n, n).T, oa_c2.reshape(n, n).T
    z_c2 = np.log10(at_c2_grid / oa_c2_grid)
    fig.add_trace(
        go.Heatmap(
            x=aas,
            y=aas,
            z=z_c1,
            colorscale="cividis",
            colorbar=dict(
                title="log<sub>10</sub>(At/Oa)",
            ),
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Heatmap(
            x=aas,
            y=aas,
            z=z_c2,
            colorscale="cividis",
            showscale=False,
        ),
        row=1,
        col=2,
    )

    fig.update_layout(
        xaxis=dict(title="Resource allocation At", dtick=0.2),
        xaxis2=dict(title="Resource allocation At", dtick=0.2),
        yaxis=dict(title="Resource allocation Oa", dtick=0.2),
        width=760,
        height=380,
        title="Coexistence in succinate-glucose space",
    )
    fig = style_plot(fig, font_size=12, line_thickness=2, top_margin=35)
    fig.write_image("plots/contours/coexistence_resource_allocation.svg")


resource_allocation_heatmap()
