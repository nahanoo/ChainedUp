from model import Model
from experiment import Species
import numpy as np
from style import *
import plotly.graph_objects as go
import pandas as pd
from os import path
from multiprocessing import Pool

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
        at.a = a[0]
        oa.a = a[1]
    model = Model(at, oa, None, xs, conc, 0.3)
    model.integrate_model()
    return model.succinate[-1], model.glucose[-1], model.at.y[-1], model.oa.y[-1]


def a_mono_culture():
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    xs = np.linspace(0, 100, 1000)
    aas = np.linspace(0.0, 1.0, 50)
    with Pool() as pool:
        args_at = [(conc, params, xs, "At", a) for a in aas]
        results_at = pool.map(simulate_endpoints, args_at)
        C1s_at, C2s_at, At_end_at, Oa_end_at = zip(*results_at)
        args_oa = [(conc, params, xs, "Oa", a) for a in aas]
        results_oa = pool.map(simulate_endpoints, args_oa)
        C1s_oa, C2s_oa, At_end_oa, Oa_end_oa = zip(*results_oa)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=C1s_at,
            y=C2s_at,
            mode="markers",
            marker=dict(
                symbol="triangle-up",
                color=aas,
                colorscale="viridis",
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
            x=C1s_oa,
            y=C2s_oa,
            mode="markers",
            marker=dict(
                color=aas,
                colorscale="viridis",
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
    fig.write_image("plots/contours/a_C1_C2_resource_space2.pdf")


def a1_a2_sweep():
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)

    xs = np.linspace(0, 100, 1000)
    aas = np.linspace(0.0, 1.0, 50)

    n = len(aas)
    Z_C1 = np.empty((n, n))
    Z_C2 = np.empty((n, n))
    Z_N1 = np.empty((n, n))
    Z_N2 = np.empty((n, n))

    for at_a in aas:
        for oa_a in aas:
            args = (conc, params, xs, "At+Oa", (at_a, oa_a))
            C1m, C2m, N1m, N2m = simulate_endpoints(args)
            j = np.where(aas == at_a)[0][0]
            i = np.where(aas == oa_a)[0][0]
            Z_C1[i, j] = C1m
            Z_C2[i, j] = C2m
            Z_N1[i, j] = N1m
            Z_N2[i, j] = N2m

    ratio = Z_C1 / (Z_C1 + Z_C2)

    Zs = [Z_C1, Z_C2, ratio]
    scale_names = ["Succinate<br>[mM]", "Glucose<br>[mM]", "Succinate<br>ratio"]
    zmaxs = [np.nanmax(Z_C1), np.nanmax(Z_C2), 1.0]  # ratio in [0,1]
    filenames = ["a1_a2_C12.pdf", "a1_a2_C22.pdf", "a1_a2_ratio2.pdf"]

    for idx, Z in enumerate(Zs):
        fig = go.Figure()
        fig.add_trace(
            go.Contour(
                z=Z,
                x=aas,
                y=aas,
                contours=dict(showlines=False),
                colorscale="Viridis",
                zmin=0,
                zmax=zmaxs[idx],
                zauto=False,
                ncontours=50,
                colorbar=dict(title=scale_names[idx], len=0.6, y=0.2, thickness=12),
                name=scale_names[idx],
            )
        )

        fig.update_layout(
            xaxis=dict(title="a At", ticks="inside"),
            yaxis=dict(title="a Oa", ticks="inside"),
        )
        fig = style_plot(fig, line_thickness=1.5, font_size=font_size)
        fig.write_image(f"plots/contours/{filenames[idx]}")

    ratio = Z_N1 / (Z_N1 + Z_N2)
    Zs = [Z_N1, Z_N2, ratio]
    scale_names = ["OD At", "OD Oa", "At ratio"]
    zmaxs = [np.nanmax(Z_N1), np.nanmax(Z_N2), 1.0]  # ratio in [0,1]
    filenames = ["a1_a2_C12_At.pdf", "a1_a2_C22_Oa.pdf", "a1_a2_ratio2_At_Oa.pdf"]

    for idx, Z in enumerate(Zs):
        fig = go.Figure()
        fig.add_trace(
            go.Contour(
                z=Z,
                x=aas,
                y=aas,
                contours=dict(showlines=False),
                colorscale="Viridis",
                zmin=0,
                zmax=zmaxs[idx],
                zauto=False,
                ncontours=50,
                colorbar=dict(title=scale_names[idx], len=0.6, y=0.2, thickness=12),
                name=scale_names[idx],
            )
        )

        fig.update_layout(
            xaxis=dict(title="a At", ticks="inside"),
            yaxis=dict(title="a Oa", ticks="inside"),
        )
        fig = style_plot(fig, line_thickness=1.5, font_size=font_size)
        fig.write_image(f"plots/contours/{filenames[idx]}")


a_mono_culture()
