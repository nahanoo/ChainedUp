import numpy as np
import pandas as pd
from os import path

import plotly.graph_objects as go
from scipy.optimize import brentq, root

from experiment import Species
from style import *


# ----------------------------
# Settings
# ----------------------------
succ_M = 3.75
gluc_M = 2.5
D = 0.3

succ_concs = np.linspace(0, succ_M, 50)
gluc_concs = np.linspace(0, gluc_M, 50)


# ----------------------------
# Core model pieces
# ----------------------------
def consumption_line(sp, succ_M, gluc_M, succ):
    """Return (gluc_vals, succ_vals) on the line, with gluc>=0 masked out."""
    m = sp.q_succ / sp.q_gluc
    succ = np.asarray(succ, dtype=float)
    gluc = gluc_M - m * (succ_M - succ)

    return gluc, succ


def mu(sp, succ, gluc):
    """Allocation-weighted Monod growth/uptake surface."""
    return sp.a * sp.v_succ * succ / (sp.K_succ + succ) + (
        1 - sp.a
    ) * sp.v_gluc * gluc / (sp.K_gluc + gluc)


def zngi_line_intersection(sp, D, n=2000):
    """Return (succ*, gluc*) where mu= D on the consumption line (gluc>=0)."""
    succ_grid = np.linspace(0, succ_M, n)
    m = sp.q_succ / sp.q_gluc
    gluc_grid = gluc_M - m * (succ_M - succ_grid)

    if succ_grid.size < 2:
        return None, None

    f = mu(sp, succ_grid, gluc_grid) - D
    idx = np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    if len(idx) == 0:
        return None, None

    i = idx[0]
    s0, s1 = succ_grid[i], succ_grid[i + 1]

    def f1(s):
        g = gluc_M - m * (succ_M - s)
        return mu(sp, s, g) - D

    s_star = brentq(f1, s0, s1)
    g_star = gluc_M - m * (succ_M - s_star)
    return float(s_star), float(g_star)


def zngi_zngi_intersection(sp1, sp2, n=150, tol=1e-4):
    """
    Return (succ*, gluc*) where mu(sp1)=D and mu(sp2)=D within the plot bounds.
    Returns (None, None) if not found.
    """
    s = np.linspace(0, succ_M, n)
    g = np.linspace(0, gluc_M, n)
    S, G = np.meshgrid(s, g)

    f1 = mu(sp1, S, G) - D
    f2 = mu(sp2, S, G) - D
    cost = f1**2 + f2**2

    j, i = np.unravel_index(np.argmin(cost), cost.shape)
    s0, g0 = float(S[j, i]), float(G[j, i])

    if float(cost[j, i]) > tol:
        return None, None

    def F(x):
        ss, gg = x
        return [mu(sp1, ss, gg) - D, mu(sp2, ss, gg) - D]

    sol = root(F, x0=[s0, g0], method="hybr")
    if not sol.success:
        return None, None

    s_star, g_star = sol.x
    if not (0 <= s_star <= succ_M and 0 <= g_star <= gluc_M):
        return None, None

    return float(s_star), float(g_star)


def coexistence_possible(sp1, sp2):
    """
    Quick graphical check (not a full dynamical proof):
      - mutual invasibility at monoculture steady states
      - ZNGI-ZNGI intersection exists
      - supply-feasible at intersection (net slope between slopes)
    """
    s1, g1 = zngi_line_intersection(sp1, D)
    s2, g2 = zngi_line_intersection(sp2, D)

    inv_2_into_1 = (s1 is not None) and (mu(sp2, s1, g1) > D)
    inv_1_into_2 = (s2 is not None) and (mu(sp1, s2, g2) > D)

    sI, gI = zngi_zngi_intersection(sp1, sp2)
    zngi_intersects = sI is not None

    feasible = False
    if zngi_intersects:
        dS = succ_M - sI
        dG = gluc_M - gI
        if abs(dS) < 1e-12:
            feasible = abs(dG) < 1e-12
        else:
            m_net = dG / dS
            m1 = sp1.q_succ / sp1.q_gluc
            m2 = sp2.q_succ / sp2.q_gluc
            feasible = min(m1, m2) <= m_net <= max(m1, m2)

    ok = inv_2_into_1 and inv_1_into_2 and zngi_intersects and feasible
    return ok


# ----------------------------
# Plotting
# ----------------------------
def plot_zngi_and_consumption(sp1, sp2, out="tmp.svg"):
    fig = go.Figure()
    Sg, Gg = np.meshgrid(succ_concs, gluc_concs)

    for sp in [sp1, sp2]:
        Z = mu(sp, Sg, Gg)
        fig.add_trace(
            go.Contour(
                x=succ_concs,
                y=gluc_concs,
                z=Z,
                contours=dict(
                    type="constraint",
                    operation="=",
                    value=D,
                    coloring="none",
                    showlines=True,
                ),
                line=dict(width=2, color=colors[sp.name]),
                showscale=False,
                name=f"{sp.name} ZNGI (μ={D})",
                hoverinfo="skip",
            )
        )

        g_line, s_line = consumption_line(sp, succ_M, gluc_M, succ_concs)
        fig.add_trace(
            go.Scatter(
                x=s_line,
                y=g_line,
                mode="lines",
                line=dict(dash="dash", color=colors[sp.name]),
                name=f"{sp.name} consumption line",
            )
        )

        s_star, g_star = zngi_line_intersection(sp, D)
        if s_star is not None:
            fig.add_trace(
                go.Scatter(
                    x=[s_star],
                    y=[g_star],
                    mode="markers",
                    marker=dict(color=colors[sp.name], size=10),
                    name=f"{sp.name} (S*,G*)",
                )
            )

    fig.add_trace(
        go.Scatter(
            x=[succ_M],
            y=[gluc_M],
            mode="markers",
            name="Supply point",
            marker=dict(size=10, color="gray"),
        )
    )

    sI, gI = zngi_zngi_intersection(sp1, sp2)
    if sI is not None:
        fig.add_trace(
            go.Scatter(
                x=[sI],
                y=[gI],
                mode="markers",
                name="ZNGI-ZNGI intersection",
                marker=dict(size=10, color="dimgray"),
            )
        )

    ok = coexistence_possible(sp1, sp2)
    fig.update_layout(
        xaxis_title="Succinate [mM]",
        yaxis_title="Glucose [mM]",
        title=f"ZNGIs + consumption lines (D={D}) | coexistence_possible={ok}",
        width=420,
        height=420,
    )

    fig = style_plot(fig, marker_size=8)
    fig.write_image(out)


def plot_invasion_map(sp1, sp2, out="invasion_map.svg"):
    a_vals = np.linspace(0.2, 0.8, 100)
    Z = np.zeros((len(a_vals), len(a_vals)), dtype=int)  # rows=y=a_sp2, cols=x=a_sp1

    # codes:
    # 0 = no invasion either way
    # 1 = sp2 invades sp1
    # 2 = sp1 invades sp2
    # 3 = mutual invasion
    for iy, a2 in enumerate(a_vals):
        for ix, a1 in enumerate(a_vals):
            sp1.a = float(a1)
            sp2.a = float(a2)

            s1, g1 = zngi_line_intersection(sp1, D)
            inv_2_into_1 = (s1 is not None) and (mu(sp2, s1, g1) > D)

            s2, g2 = zngi_line_intersection(sp2, D)
            inv_1_into_2 = (s2 is not None) and (mu(sp1, s2, g2) > D)

            Z[iy, ix] = (1 if inv_2_into_1 else 0) + (2 if inv_1_into_2 else 0)

    labels = {
        0: "no invasion",
        1: f"{sp2.name} invades {sp1.name}",
        2: f"{sp1.name} invades {sp2.name}",
        3: "mutual invasion",
    }

    colorscale = [
        [0.00, "#f2f2f2"],
        [0.24, "#f2f2f2"],  # 0
        [0.25, "#99bbff"],
        [0.49, "#99bbff"],  # 1
        [0.50, "#ff9999"],
        [0.74, "#ff9999"],  # 2
        [0.75, "#b6f2c2"],
        [1.00, "#b6f2c2"],  # 3
    ]

    fig = go.Figure(
        go.Heatmap(
            x=a_vals,
            y=a_vals,
            z=Z,
            zmin=0,
            zmax=3,
            colorscale=colorscale,
            colorbar=dict(
                tickmode="array",
                tickvals=[0, 1, 2, 3],
                ticktext=[labels[k] for k in [0, 1, 2, 3]],
            ),
            hovertemplate=f"a_{sp1.name}=%{{x:.2f}}<br>a_{sp2.name}=%{{y:.2f}}<br>code=%{{z}}<extra></extra>",
        )
    )

    fig.update_layout(
        xaxis_title=f"a ({sp1.name})",
        yaxis_title=f"a ({sp2.name})",
        title=f"Invasion map at D={D}",
        width=520,
        height=480,
    )

    fig = style_plot(fig, marker_size=8)
    fig.write_image(out)


# ----------------------------
# Run
# ----------------------------
conc_C = 15
p_f = path.join("parameters", f"parameters_{conc_C}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)

at = Species("At", params.loc["At"])
oa = Species("Oa", params.loc["Oa"])

at.a = 0.1
oa.a = 0.7

plot_zngi_and_consumption(at, oa, out="tmp.svg")
plot_invasion_map(at, oa, out="invasion_map.svg")
