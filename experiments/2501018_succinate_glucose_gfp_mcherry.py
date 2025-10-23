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


def strip_numbers(y, threshold=100):
    stripped = []
    for i in y:
        s = str(int(i)).rstrip("0")
        if s == "":  # happens for values like 0, 00, 000
            s = "0"
        stripped.append(int(s))
    return stripped


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


species = ["At", "Oa", "At+Oa"]
c_mono = 30
concs = [c_mono, c_mono, 2 * c_mono]
css = ["Succinate", "Glucose", "Succinate+Glucose"]
titles = ["OD600", "GFP", "mCherry"]
signal_OD = [
    "2510_succinate_glucose_251018_succinate",
    "2510_succinate_glucose_251018_glucose",
    "2510_succinate_glucose_251018_succinate_glucose",
]
signals_gfp = [
    "2510_succinate_glucose_251018_succinate_gfp",
    "2510_succinate_glucose_251018_glucose_gfp",
    "2510_succinate_glucose_251018_succinate_glucose_gfp",
]
signals_mcherry = [
    "2510_succinate_glucose_251018_succinate_mcherry",
    "2510_succinate_glucose_251018_glucose_mcherry",
    "2510_succinate_glucose_251018_succinate_glucose_mcherry",
]
fs = [
    "plots/growth_curves/succinate_gfp_mcherry.pdf",
    "plots/growth_curves/glucose_gfp_mcherry.pdf",
    "plots/growth_curves/succinate_glucose_gfp_mcherry.pdf",
]

colors = viridis_palette_plotly(titles)
for i, cs in enumerate(css):
    fig = make_subplots(rows=3, cols=1, subplot_titles=titles)
    for j, signal in enumerate([signal_OD[i], signals_gfp[i], signals_mcherry[i]]):
        linegroups = get_linegroups("At+Oa", cs, concs[i], meta, signal)
        for lg in linegroups:
            x = data[lg + "_time"].to_numpy()
            y = data[lg + "_measurement"].to_numpy()
            x, y = filter_time(x, y, 35)
            if signal in signals_gfp + signals_mcherry:
                y = strip_numbers(y)
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    line=dict(color=colors[titles[i]]),
                    name=titles[i],
                    showlegend=False,
                ),
                col=1,
                row=j + 1,
            )
    fig.update_layout(
        title="Co-culture on " + cs,
        xaxis3_title="Time (hours)",
        yaxis_title="OD600",
        yaxis2_title="GFP",
        yaxis3_title="mCherry",
    )
    fig = style_plot(fig, font_size=font_size)
    fig.write_image(fs[i])
