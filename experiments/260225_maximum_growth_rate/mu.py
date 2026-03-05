import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from parse_data import parse_at, parse_oa
from style import *  # expects: colors, style_plot


# ----------------------------
# Config
# ----------------------------
SPECIES = ["At", "Oa"]
SUBSTRATES = ["Glucose", "Succinate"]

SUBPLOT_TITLES = {
    ("At", "Glucose"): "At Glucose",
    ("At", "Succinate"): "At Succinate",
    ("Oa", "Glucose"): "Oa Glucose",
    ("Oa", "Succinate"): "Oa Succinate",
}

OUT_GROWTH = "growth_curves.svg"
OUT_MU = "max_growth_rates.svg"

WINDOW_H = 3.0  # sliding window size in hours
PEAK_MODE = "max"  # "max" or "topk_mean"
TOPK = 5  # used if PEAK_MODE == "topk_mean"


# ----------------------------
# Helpers
# ----------------------------
def make_2x2():
    return make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            SUBPLOT_TITLES[(s, sub)] for s in SPECIES for sub in SUBSTRATES
        ],
        vertical_spacing=0.10,
        horizontal_spacing=0.10,
    )


def panel_rc(i: int):
    """Map i=0..3 to (row, col) in a 2x2 grid."""
    return (i // 2 + 1, i % 2 + 1)


def subset_condition(df: pd.DataFrame, sp: str, sub: str) -> pd.DataFrame:
    return df[(df["species"] == sp) & (df["substrate"] == sub)].copy()


def compute_mu_sliding(rep: pd.DataFrame, window_h: float) -> pd.DataFrame:
    """
    Sliding-window estimate of instantaneous growth rate:
      mu = slope of ln(OD) ~ t on [t_i, t_i + window_h]
    Returns DF with columns: timepoint, mu
    """
    rep = rep.sort_values("timepoint")

    x = rep["timepoint"].to_numpy(dtype=float)
    y = rep["OD"].to_numpy(dtype=float)

    # log requires positive OD
    ok = np.isfinite(x) & np.isfinite(y) & (y > 0)
    x, y = x[ok], y[ok]

    mus, t_mids = [], []
    for t0 in x:
        t1 = t0 + window_h
        m = (x >= t0) & (x <= t1)
        if m.sum() < 2:
            continue

        xw = x[m]
        yw = y[m]
        if np.any(yw <= 0):
            continue

        mu = np.polyfit(xw, np.log(yw), 1)[0]
        mus.append(mu)
        t_mids.append(xw.mean())

    return pd.DataFrame({"timepoint": t_mids, "mu": mus})


def peak_mu(mus: np.ndarray, mode: str = "max", topk: int = 5) -> float:
    mus = np.asarray(mus, dtype=float)
    mus = mus[np.isfinite(mus)]
    if len(mus) == 0:
        return np.nan

    if mode == "max":
        return float(np.max(mus))

    if mode == "topk_mean":
        k = min(topk, len(mus))
        return float(np.mean(np.sort(mus)[-k:]))

    raise ValueError(f"Unknown mode: {mode}")


# ----------------------------
# Load data once
# ----------------------------
df = pd.concat([parse_at(), parse_oa()], ignore_index=True)
df = df[["species", "substrate", "replicate", "timepoint", "OD"]].copy()
df = df.dropna(subset=["species", "substrate", "replicate", "timepoint", "OD"])


# ----------------------------
# 1) Growth curves
# ----------------------------
fig = make_2x2()

for k, (sp, sub) in enumerate([(s, c) for s in SPECIES for c in SUBSTRATES]):
    row, col = panel_rc(k)
    subdf = subset_condition(df, sp, sub)

    for rep_id, rep in subdf.groupby("replicate", sort=True):
        rep = rep.sort_values("timepoint")

        fig.add_trace(
            go.Scatter(
                x=rep["timepoint"],
                y=rep["OD"],
                mode="markers",
                marker=dict(
                    color=colors[sp],
                    opacity=0.55,
                ),
                name=f"{sp} {sub} {rep_id}",
                showlegend=False,
            ),
            row=row,
            col=col,
        )

# axis labels
for c in [1, 2]:
    fig.update_xaxes(title_text="Time (h)", row=2, col=c)
for r in [1, 2]:
    fig.update_yaxes(title_text="OD", row=r, col=1)

fig = style_plot(fig, font_size=12, marker_size=4)
fig.write_image(OUT_GROWTH)


# ----------------------------
# 2) Sliding-window growth rates (mu)
# ----------------------------
fig = make_2x2()

for k, (sp, sub) in enumerate([(s, c) for s in SPECIES for c in SUBSTRATES]):
    row, col = panel_rc(k)
    subdf = subset_condition(df, sp, sub)

    rep_peaks = []

    for rep_id, rep in subdf.groupby("replicate", sort=True):
        mu_df = compute_mu_sliding(rep, window_h=WINDOW_H)
        if len(mu_df) == 0:
            continue

        # store peak per replicate (not mean!)
        rep_peaks.append(peak_mu(mu_df["mu"].to_numpy(), mode=PEAK_MODE, topk=TOPK))

        fig.add_trace(
            go.Scatter(
                x=mu_df["timepoint"],
                y=mu_df["mu"],
                mode="markers",
                marker=dict(color=colors[sp], opacity=0.55),
                name=f"{sp} {sub} {rep_id}",
                showlegend=False,
            ),
            row=row,
            col=col,
        )

    # one annotation per panel
    if len(rep_peaks) > 0 and np.isfinite(rep_peaks).any():
        mu_max = float(np.nanmax(rep_peaks))
        fig.add_annotation(
            text=f"Max μ ({PEAK_MODE}): {mu_max:.3f} h⁻¹",
            x=0.98,
            y=0.02,
            xref="paper",
            yref="paper",
            xanchor="right",
            yanchor="bottom",
            showarrow=False,
            row=row,
            col=col,
        )

# axis labels
for c in [1, 2]:
    fig.update_xaxes(title_text="Time (h)", row=2, col=c)
for r in [1, 2]:
    fig.update_yaxes(title_text="μ (h⁻¹)", row=r, col=1)

fig = style_plot(fig, font_size=12, marker_size=4)
fig.write_image(OUT_MU)
