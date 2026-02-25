import os
import os.path as path
from experiment import Species
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from style import *  # style_plot, colors, etc.

# from your_module import Species  # assumes you already have this


# ---------------------------
# Model helpers
# ---------------------------
C_SUCC = 4  # carbons in succinate
C_GLUC = 6  # carbons in glucose


def growth_rate(sp, s, g):
    """Your J(s,g): allocation-weighted Monod uptake/growth surface."""
    return sp.a * sp.v_succ * s / (sp.K_succ + s) + (1 - sp.a) * sp.v_gluc * g / (
        sp.K_gluc + g
    )


def get_q_gluc(sp):
    """Support both q_glucc (your attribute) and q_gluc (common typo)."""
    if hasattr(sp, "q_glucc"):
        return sp.q_glucc
    return sp.q_gluc  # will raise if neither exists (good: fail loudly)


def chemostat_inflow_from_C(conc_C, frac_succ_C):
    """
    Build inflow point (S_in, G_in) from total carbon concentration conc_C [mM C]
    and fraction of carbon coming from succinate.
    """
    S_in = (frac_succ_C * conc_C) / C_SUCC
    G_in = ((1 - frac_succ_C) * conc_C) / C_GLUC
    return S_in, G_in


def consumption_line(sp, S_in, G_in, s_vals):
    """
    Supply/consumption line through (S_in, G_in) with slope m = q_succ / q_gluc.
    G = G_in - m*(S_in - S)
    """
    m = sp.q_succ / get_q_gluc(sp)
    g_vals = G_in - m * (S_in - s_vals)
    return s_vals, g_vals


def predict_resource_ss(sp, D, S_in, G_in, n=2000):
    """
    Predict (S*,G*) by intersecting the line with mu(s,g)=D.
    We trace the line from S=0..S_in, compute mu, and find the first sign change of mu-D.
    """
    s = np.linspace(0, S_in, n)
    s, g = consumption_line(sp, S_in, G_in, s)

    # keep physical region
    msk = (g >= 0) & np.isfinite(g)
    s, g = s[msk], g[msk]
    if len(s) < 2:
        return None, None

    mu = growth_rate(sp, s, g)
    f = mu - D

    # find bracketing interval for root
    sign = np.sign(f)
    idx = np.where(sign[:-1] * sign[1:] <= 0)[0]
    if len(idx) == 0:
        return None, None

    i = idx[0]
    # linear interpolation
    s0, s1 = s[i], s[i + 1]
    g0, g1 = g[i], g[i + 1]
    f0, f1 = f[i], f[i + 1]

    if f1 == f0:
        return float(s0), float(g0)

    w = -f0 / (f1 - f0)
    s_star = s0 + w * (s1 - s0)
    g_star = g0 + w * (g1 - g0)
    return float(s_star), float(g_star)


# ---------------------------
# Load parameters + species
# ---------------------------
conc_C = 15
p_f = path.join("parameters", f"parameters_{conc_C}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)

at = Species("At", params.loc["At"])
oa = Species("Oa", params.loc["Oa"])

# allocation parameters
at.a = 0.2
oa.a = 0.8

# ---------------------------
# Chemostat settings
# ---------------------------
D = 0.3  # h^-1, the isocline level you are plotting

# choose inflow point (S_in, G_in)
# Option A: define by carbon split of total conc_C (recommended)
frac_succ_C = 0.5  # 50% of carbon from succinate, 50% from glucose
S_in, G_in = chemostat_inflow_from_C(conc_C * 2, frac_succ_C)

# Option B: or set directly (uncomment)
# S_in, G_in = 1.5, 0.8

# plotting ranges consistent with your conc_C definition
s_concs = np.linspace(0, conc_C / C_SUCC, 250)  # 0..3.75
g_concs = np.linspace(0, conc_C / C_GLUC, 250)  # 0..2.5
s_grid, g_grid = np.meshgrid(s_concs, g_concs)

JAt_grid = growth_rate(at, s_grid, g_grid)
JOa_grid = growth_rate(oa, s_grid, g_grid)

# ---------------------------
# Figure: ZNGIs + consumption lines + predicted (S*,G*)
# ---------------------------
fig = go.Figure()


def add_zngi(fig, name, zgrid, color):
    fig.add_trace(
        go.Contour(
            x=s_concs,
            y=g_concs,
            z=zgrid,
            contours=dict(
                type="constraint",
                operation="=",
                value=D,
                coloring="none",
                showlines=True,
            ),
            line=dict(width=2, color=color),
            name=f"{name}: μ=D",
            showscale=False,
            hoverinfo="skip",
            showlegend=True,
        )
    )


add_zngi(fig, "At", JAt_grid, "red")
add_zngi(fig, "Oa", JOa_grid, "blue")

# supply point
fig.add_trace(
    go.Scatter(
        x=[S_in],
        y=[G_in],
        mode="markers",
        marker=dict(size=10),
        name="Inflow (S_in, G_in)",
    )
)


def add_consumption_and_ss(fig, sp, color):
    s_line = np.linspace(0, max(s_concs), 400)
    s_line, g_line = consumption_line(sp, S_in, G_in, s_line)
    msk = (g_line >= 0) & (g_line <= max(g_concs))
    fig.add_trace(
        go.Scatter(
            x=s_line[msk],
            y=g_line[msk],
            mode="lines",
            line=dict(dash="dash", width=2, color=color),
            name=f"{sp.name}: consumption line",
        )
    )

    s_star, g_star = predict_resource_ss(sp, D, S_in, G_in)
    if (s_star is not None) and (g_star is not None):
        fig.add_trace(
            go.Scatter(
                x=[s_star],
                y=[g_star],
                mode="markers+text",
                marker=dict(size=9, color=color),
                text=[f"{sp.name} (S*,G*)"],
                textposition="top right",
                name=f"{sp.name}: predicted (S*,G*)",
            )
        )


add_consumption_and_ss(fig, at, "red")
add_consumption_and_ss(fig, oa, "blue")

fig.update_layout(
    xaxis_title="Succinate [mM]",
    yaxis_title="Glucose [mM]",
    title=f"Dual-substrate chemostat: μ=D (D={D}) + inflow & predicted resource steady state",
    width=420,
    height=420,
)

fig = style_plot(fig, font_size=14, line_thickness=2)

os.makedirs(path.join("plots", "contours"), exist_ok=True)
fig.write_image("plots/contours/isoclines.svg")
