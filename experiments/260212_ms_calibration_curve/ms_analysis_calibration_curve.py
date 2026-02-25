import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.colors import qualitative
from style import *

path = "../../data/260212_ms_calibration_curve/Quant_results_area.csv"
raw = pd.read_csv(path)

df = raw.iloc[1:].copy()
df = df.rename(
    columns={
        "Unnamed: 2": "Name",
        "Unnamed: 4": "Type",
        "Unnamed: 13": "Sample_type",
        "Unnamed: 14": "Conc_mM",
        "Unnamed: 1": "Flag",
        "Unnamed: 15": "Well",
    }
)

num_cols = [
    "Conc_mM",
    "Norvaline Results",
    "Norleucine Results",
    "Sorbitol Results",
    "Succinate Results",
    "Glucose_1 Results",
    "Glucose_2 Results",
]
for c in num_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df["Glucose_area"] = df["Glucose_1 Results"] + df["Glucose_2 Results"]

df["Glucose_resp"] = df["Glucose_area"] / df["Norvaline Results"]
df["Succinate_resp"] = df["Succinate Results"] / df["Sorbitol Results"]


glc = df[(df["Type"] == "Sample") & (df["Sample_type"] == "Glucose")].copy()
suc = df[(df["Type"] == "Sample") & (df["Sample_type"] == "Succinate")].copy()
glc_concs = glc["Conc_mM"].sort_values()
suc_concs = suc["Conc_mM"].sort_values()
glc = glc.sort_values("Conc_mM")
suc = suc.sort_values("Conc_mM")
blank = df[(df["Type"] == "Blank") & (df["Sample_type"] == "Blank")].copy()


def lin_reg(x, y):
    coeffs = np.polyfit(x, y, 1, w=1 / np.sqrt(x))
    return coeffs


glc_coeffs = lin_reg(glc["Conc_mM"], glc["Glucose_resp"])
suc_coeffs = lin_reg(suc["Conc_mM"], suc["Succinate_resp"])
fit_glc = np.poly1d(glc_coeffs)
fit_suc = np.poly1d(suc_coeffs)
x_g = np.logspace(np.log10(glc["Conc_mM"].min()), np.log10(glc["Conc_mM"].max()), 100)
x_s = np.logspace(np.log10(suc["Conc_mM"].min()), np.log10(suc["Conc_mM"].max()), 100)

colors = {
    "Glucose": qualitative.Bold[0],
    "Succinate": qualitative.Bold[1],
    "Blank Succinate": qualitative.Bold[2],
    "Blank Glucose": qualitative.Bold[3],
}


def get_concentration(response, fit):
    return (response - fit.coeffs[1]) / fit.coeffs[0]


glc_conc = [get_concentration(r, fit_glc) for r in glc["Glucose_resp"]]
suc_conc = [get_concentration(r, fit_suc) for r in suc["Succinate_resp"]]


fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=glc["Conc_mM"],
        y=glc["Glucose_resp"],
        mode="markers",
        textposition="top center",
        name="Glucose",
        text=glc["Well"],
        marker=dict(color=colors["Glucose"], line=dict(color="black", width=1.2)),
    )
)
"""fig.add_trace(
    go.Scatter(
        x=x_g,
        y=fit_glc(x_g),
        mode="lines",
        name="Glucose Fit",
        line=dict(color=colors["Glucose"]),
    )
)"""
fig.add_trace(
    go.Scatter(
        x=suc["Conc_mM"],
        y=suc["Succinate_resp"],
        mode="markers",
        text=suc["Well"],
        textposition="top center",
        name="Succinate",
        marker=dict(color=colors["Succinate"], line=dict(color="black", width=1.2)),
    )
)
"""fig.add_trace(
    go.Scatter(
        x=x_s,
        y=fit_suc(x_s),
        mode="lines",
        name="Succinate Fit",
        line=dict(color=colors["Succinate"]),
    )
)"""

fig.update_xaxes(
    title="Concentration (mM)",
    ticks="inside",
    type="log",
    exponentformat="power",
)
fig.update_yaxes(
    title="Area",
    ticks="inside",
    exponentformat="power",
    type="log",
)
fig.update_layout(
    width=600,
    height=400,
    legend=dict(
        yanchor="top",
        y=0.95,
        xanchor="left",
        x=0.05,
        bgcolor="rgba(255,255,255,1)",
    ),
    title="Calibration Curves for Glucose and Succinate",
)
fig = style_plot(fig, marker_size=10, line_thickness=2, font_size=14)
fig.write_image("plots/ms_calibration_curve.pdf")
