import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from plotly.subplots import make_subplots
import plotly.express as px


meta = pd.read_csv("../data/251018_succinate_glucose_plate_reader/metaod/metadata.csv")
data = pd.read_csv(
    "../data/251018_succinate_glucose_plate_reader/metaod/measurements.csv",
    index_col=0,
)
colors = {"Succinate": "#31688e", "Glucose": "#35b779", "Succinate+Glucose": "#fde725"}

font_size = 14


def filter_time(x, y, cut_off):
    mask = x <= cut_off
    return x[mask], y[mask]


def viridis_palette_plotly(cs):
    xs = np.linspace(0.12, 0.88, len(cs))
    cols = px.colors.sample_colorscale("Viridis", xs)
    return dict(zip(cs, cols))


def get_linegroups(species, cs, conc, meta, signal):
    linegroups = meta[
        (meta["species"] == species)
        & (meta["carbon_source"] == cs)
        & (meta["cs_conc"] == conc)
        & (meta["exp_ID"] == signal)
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
    # build a color map keyed by concentration values (e.g., 1, 2, 5 mM)
    concs = list(dict.fromkeys(C_to_mM_succinate.values()))  # preserves order
    colors = viridis_palette_plotly(concs)

    fig = make_subplots(
        rows=3, cols=1, subplot_titles=species, shared_xaxes=True, vertical_spacing=0.08
    )
    seen_conc = set()  # which concentrations already have a legend entry?

    for s in species:
        for c, conc in C_to_mM_succinate.items():
            linegroups = get_linegroups(
                s, "Succinate", c, meta, "2510_succinate_glucose_251018_succinate"
            )
            for lg in linegroups:
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()

                x, y = filter_time(x, y, 35)  # Apply time filtering

                show = conc not in seen_conc  # <- only first of this conc shows
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        name=f"{conc} mM<br>Succinate",
                        legendgroup=str(conc),  # group by concentration
                        showlegend=show,
                        line=dict(color=colors[conc], width=2),
                        # legendgrouptitle_text=f"{conc} mM"  # optional: group title
                    ),
                    row=species.index(s) + 1,
                    col=1,
                )
                if show:
                    seen_conc.add(conc)

    fig.update_layout(
        title="Growth on Succinate",
        xaxis3_title="Time (hours)",
        yaxis2_title="OD600",
    )
    fig = style_plot(fig, font_size=font_size)
    fig.write_image("plots/growth_curves/succinate.pdf")


def glucose():
    # colors keyed by glucose concentration (e.g., 0.5, 1, 2 mM)
    concs = list(dict.fromkeys(C_to_mM_glucose.values()))
    colors = viridis_palette_plotly(concs)

    fig = make_subplots(
        rows=3, cols=1, subplot_titles=species, shared_xaxes=True, vertical_spacing=0.08
    )
    seen_conc = set()

    for s in species:
        for c, conc in C_to_mM_glucose.items():
            linegroups = get_linegroups(
                s, "Glucose", c, meta, "2510_succinate_glucose_251018_glucose"
            )
            for lg in linegroups:
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                x, y = filter_time(x, y, 35)  # Apply time filtering
                show = conc not in seen_conc
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        name=f"{conc} mM<br>Glucose",
                        legendgroup=str(conc),
                        showlegend=show,
                        line=dict(color=colors[conc], width=2),
                    ),
                    row=species.index(s) + 1,
                    col=1,
                )
                if show:
                    seen_conc.add(conc)

    fig.update_layout(
        title="Growth on Glucose",
        xaxis3_title="Time (hours)",
        yaxis2_title="OD600",
    )
    fig = style_plot(fig, font_size=font_size)
    fig.write_image("plots/growth_curves/glucose.pdf")


def glucose_succinate():
    """
    Legend shows each unique (Succinate, Glucose) pair once.
    Colors keyed by the pair (as tuples).
    """
    # unique (succ, glc) pairs as tuples, preserving order
    pairs_all = [tuple(v) for v in C_to_mM_glucose_succinate.values()]
    pairs = list(dict.fromkeys(pairs_all))  # dedupe while keeping order
    colors = viridis_palette_plotly(pairs)  # keys are tuples

    fig = make_subplots(
        rows=3, cols=1, subplot_titles=species, shared_xaxes=True, vertical_spacing=0.08
    )
    seen_pair = set()

    for s in species:
        for c, conc_pair in C_to_mM_glucose_succinate.items():
            succ_mM, glc_mM = tuple(conc_pair)  # make sure it's a tuple
            linegroups = get_linegroups(
                s,
                "Succinate+Glucose",
                c,
                meta,
                "2510_succinate_glucose_251018_succinate_glucose",
            )
            for lg in linegroups:
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                x, y = filter_time(x, y, 35)  # Apply time filtering
                pair_key = (succ_mM, glc_mM)
                label = f"{succ_mM} mM Succinate<br> {glc_mM} mM Glucose"
                show = pair_key not in seen_pair

                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode="lines",
                        name=label,
                        legendgroup=str(pair_key),
                        showlegend=show,
                        line=dict(color=colors[pair_key], width=2),
                    ),
                    row=species.index(s) + 1,
                    col=1,
                )
                if show:
                    seen_pair.add(pair_key)

    fig.update_layout(
        title="Growth on Succinate+Glucose",
        xaxis3_title="Time (hours)",
        yaxis2_title="OD600",
    )
    fig = style_plot(fig, font_size=font_size)
    fig.write_image("plots/growth_curves/glucose_succinate.pdf")


def glucose_succinate_layered():
    fig = make_subplots(
        rows=3, cols=1, subplot_titles=species, shared_xaxes=True, vertical_spacing=0.08
    )
    css = ["Glucose", "Succinate", "Succinate+Glucose"]
    signals = [
        "2510_succinate_glucose_251018_glucose",
        "2510_succinate_glucose_251018_succinate",
        "2510_succinate_glucose_251018_succinate_glucose",
    ]
    q = 30
    concs = [q, q, 2 * q]
    seen = set()
    title = (
        "Growth on "
        + str(C_to_mM_glucose_succinate[q * 2][1])
        + " mM Glucose and "
        + str(C_to_mM_glucose_succinate[q * 2][0])
        + " mM Succinate"
    )
    for i, cs in enumerate(css):
        for j, s in enumerate(species):
            linegroups = get_linegroups(s, cs, concs[i], meta, signals[i])
            for lg in linegroups:
                x = data[lg + "_time"].to_numpy()
                y = data[lg + "_measurement"].to_numpy()
                x, y = filter_time(x, y, 35)  # Apply time filtering
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
                    row=j + 1,
                    col=1,
                )
                if show:
                    seen.add(cs)
    fig.update_layout(
        title=title,
        xaxis3_title="Time (hours)",
        yaxis2_title="OD600",
    )
    fig = style_plot(fig, font_size=font_size)
    fig.write_image(
        "plots/growth_curves/glucose_succinate_" + title.replace(" ", "_") + ".pdf"
    )


# glucose()
# succinate()
# glucose_succinate()
glucose_succinate_layered()
