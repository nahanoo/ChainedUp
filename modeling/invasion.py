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

# IMPORTANT:
# This must match how your chemostat simulation consumes resources.
# - "split"      : dS uptake uses mu_s/q_succ and dG uses mu_g/q_gluc  (recommended for your additive model)
# - "total"      : both resources consumed proportional to mu_total (often accidentally used)
# - "alloc_total": consumed proportional to a*mu_total and (1-a)*mu_total
UPTAKE_MODE = "split"

A_EPS = 1e-10


# ----------------------------
# Core model pieces
# ----------------------------
def mu_parts(sp, succ, gluc):
    """Return (mu_succ, mu_gluc, mu_total)."""
    mu_s = sp.a * sp.v_succ * succ / (sp.K_succ + succ)
    mu_g = (1 - sp.a) * sp.v_gluc * gluc / (sp.K_gluc + gluc)
    return mu_s, mu_g, mu_s + mu_g


def mu(sp, succ, gluc):
    """Allocation-weighted Monod growth surface."""
    return mu_parts(sp, succ, gluc)[2]


def uptake_rates(sp, succ, gluc, mode=UPTAKE_MODE):
    """
    Return (u_succ, u_gluc) = resource uptake rates per biomass [mM / (OD*h)].
    This is the ONLY place you should encode how the chemostat consumes resources.
    """
    mu_s, mu_g, mu_t = mu_parts(sp, succ, gluc)

    if mode == "split":
        # each resource pays only for its own growth contribution
        u_s = mu_s / sp.q_succ
        u_g = mu_g / sp.q_gluc
    elif mode == "total":
        # both resources pay for total growth (keeps straight lines)
        u_s = mu_t / sp.q_succ
        u_g = mu_t / sp.q_gluc
    elif mode == "alloc_total":
        # allocation scales which resource pays for total growth
        u_s = sp.a * mu_t / sp.q_succ
        u_g = (1 - sp.a) * mu_t / sp.q_gluc
    else:
        raise ValueError("UPTAKE_MODE must be 'split', 'total', or 'alloc_total'")

    return u_s, u_g


def balance_residual(sp, succ, gluc, mode=UPTAKE_MODE):
    """
    Monoculture chemostat steady-state consistency (eliminating biomass X):

        D(S_in - S) = u_s(S,G) * X
        D(G_in - G) = u_g(S,G) * X

    => (S_in - S) * u_g(S,G) - (G_in - G) * u_s(S,G) = 0
    """
    u_s, u_g = uptake_rates(sp, succ, gluc, mode=mode)
    return (succ_M - succ) * u_g - (gluc_M - gluc) * u_s


# ----------------------------
# "Consumption line" (actually a curve for additive uptake)
# ----------------------------
def consumption_line(sp, succ_vals, mode=UPTAKE_MODE):
    """
    For additive uptake, the "consumption line" is a CURVE defined by:
        balance_residual(sp, S, G) = 0
    Parameterize by S and solve for G in [0, gluc_M] with brentq (when possible).

    Returns (gluc_curve, succ_curve) arrays.
    """
    succ_vals = np.asarray(succ_vals, dtype=float)

    # succ-only: glucose stays at inflow
    if sp.a >= 1.0 - A_EPS:
        return np.full_like(succ_vals, gluc_M), succ_vals

    # gluc-only: succinate stays at inflow (plot a vertical line)
    if sp.a <= A_EPS:
        g = np.linspace(gluc_M, 0.0, succ_vals.size)
        s = np.full_like(g, succ_M)
        return g, s

    gluc_curve = []
    succ_curve = []

    for s in succ_vals:

        def f(g):
            return balance_residual(sp, s, g, mode=mode)

        f0 = f(0.0)
        f1 = f(gluc_M)

        # boundary roots
        if np.isfinite(f0) and abs(f0) < 1e-12:
            g = 0.0
        elif np.isfinite(f1) and abs(f1) < 1e-12:
            g = gluc_M
        elif np.sign(f0) * np.sign(f1) < 0:
            g = brentq(f, 0.0, gluc_M)
        else:
            continue

        if 0.0 <= g <= gluc_M:
            gluc_curve.append(g)
            succ_curve.append(s)

    return np.asarray(gluc_curve, dtype=float), np.asarray(succ_curve, dtype=float)


# ----------------------------
# Resident monoculture steady state for invasion tests
# ----------------------------
def resident_ss(sp, D, n_curve=300, tol_hit=1e-5, mode=UPTAKE_MODE):
    """
    Return (succ*, gluc*, alive, reason)

    - alive=True  means non-washout (X*>0) steady state found in bounds.
    - alive=False means washout steady state (X*=0), which ALWAYS exists:
        (succ*, gluc*) = (succ_M, gluc_M)

    We find the non-washout SS by intersecting:
      mu(S,G)=D   AND   balance_residual(S,G)=0

    If numerical bracketing misses (tangency), we accept the closest point
    on the curve if abs(mu-D) is small.
    """
    mu_in = mu(sp, succ_M, gluc_M)
    if mu_in <= D:
        return float(succ_M), float(gluc_M), False, "washout"

    # a≈1: effectively succ-only -> analytic S*
    if sp.a >= 1.0 - A_EPS:
        if sp.v_succ <= D:
            return float(succ_M), float(gluc_M), False, "washout_succ_only"
        s_star = D * sp.K_succ / (sp.v_succ - D)
        if 0 <= s_star <= succ_M:
            return float(s_star), float(gluc_M), True, "succ_only"
        return float(succ_M), float(gluc_M), False, "washout_succ_only_outside"

    # a≈0: effectively gluc-only -> analytic G*
    if sp.a <= A_EPS:
        if sp.v_gluc <= D:
            return float(succ_M), float(gluc_M), False, "washout_gluc_only"
        g_star = D * sp.K_gluc / (sp.v_gluc - D)
        if 0 <= g_star <= gluc_M:
            return float(succ_M), float(g_star), True, "gluc_only"
        return float(succ_M), float(gluc_M), False, "washout_gluc_only_outside"

    # mass-balance curve
    succ_grid = np.linspace(0, succ_M, n_curve)
    g_curve, s_curve = consumption_line(sp, succ_grid, mode=mode)
    if s_curve.size < 2:
        return float(succ_M), float(gluc_M), False, "curve_missing"

    # sort for interpolation
    order = np.argsort(s_curve)
    s_curve = s_curve[order]
    g_curve = g_curve[order]

    f = mu(sp, s_curve, g_curve) - D

    # direct hit
    j = int(np.argmin(np.abs(f)))
    if abs(f[j]) <= tol_hit:
        return float(s_curve[j]), float(g_curve[j]), True, "curve_closest"

    # sign change
    idx = np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    if len(idx) == 0:
        # no sign change -> likely tangency or mismatch; return washout to avoid "no SS"
        return float(succ_M), float(gluc_M), False, "no_bracket->washout"

    i = int(idx[0])
    s0, s1 = float(s_curve[i]), float(s_curve[i + 1])

    def f1(s):
        g = float(np.interp(s, s_curve, g_curve))
        return mu(sp, s, g) - D

    s_star = brentq(f1, s0, s1)
    g_star = float(np.interp(s_star, s_curve, g_curve))
    return float(s_star), float(g_star), True, "brentq"


# ----------------------------
# ZNGI-ZNGI intersection (for plotting)
# ----------------------------
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
    Quick check based on mutual invasibility using RESIDENT steady states:
      - Oa invades At if mu(Oa at At-SS) > D
      - At invades Oa if mu(At at Oa-SS) > D
    """
    s1, g1, alive1, _ = resident_ss(sp1, D)
    s2, g2, alive2, _ = resident_ss(sp2, D)

    if not alive1 or not alive2:
        return False

    inv_2_into_1 = mu(sp2, s1, g1) > D
    inv_1_into_2 = mu(sp1, s2, g2) > D
    return inv_2_into_1 and inv_1_into_2


# ----------------------------
# Plotting
# ----------------------------
def plot_zngi_and_consumption(sp1, sp2, out="tmp.svg"):
    fig = go.Figure()
    Sg, Gg = np.meshgrid(succ_concs, gluc_concs)

    for sp in [sp1, sp2]:
        # ZNGI: mu = D
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

        # mass-balance curve ("consumption line")
        g_line, s_line = consumption_line(
            sp, np.linspace(0, succ_M, 300), mode=UPTAKE_MODE
        )
        fig.add_trace(
            go.Scatter(
                x=s_line,
                y=g_line,
                mode="lines",
                line=dict(dash="dash", color=colors[sp.name]),
                name=f"{sp.name} mass-balance curve",
            )
        )

        # resident SS marker (or washout)
        s_star, g_star, alive, reason = resident_ss(sp, D, mode=UPTAKE_MODE)
        fig.add_trace(
            go.Scatter(
                x=[s_star],
                y=[g_star],
                mode="markers",
                marker=dict(
                    color=colors[sp.name], size=10, symbol="x" if alive else "circle"
                ),
                name=f"{sp.name} SS ({'alive' if alive else 'washout'})",
                hovertemplate=f"{sp.name}<br>S*=%{{x:.3f}}<br>G*=%{{y:.3f}}<br>{reason}<extra></extra>",
            )
        )

    # supply point
    fig.add_trace(
        go.Scatter(
            x=[succ_M],
            y=[gluc_M],
            mode="markers",
            name="Supply point",
            marker=dict(size=10, color="gray"),
        )
    )

    # ZNGI-ZNGI intersection marker
    sI, gI = zngi_zngi_intersection(sp1, sp2)
    if sI is not None:
        fig.add_trace(
            go.Scatter(
                x=[sI],
                y=[gI],
                mode="markers",
                name="ZNGI-ZNGI intersection",
                marker=dict(size=10, color="dimgray", symbol="diamond"),
            )
        )

    ok = coexistence_possible(sp1, sp2)
    fig.update_layout(
        xaxis_title="Succinate [mM]",
        yaxis_title="Glucose [mM]",
        title=f"ZNGIs + mass-balance curves (mode={UPTAKE_MODE}) | coexistence_possible={ok}",
        width=420,
        height=420,
    )

    fig = style_plot(fig, marker_size=8)
    fig.write_image(out)


def plot_invasion_map(sp1, sp2, out="invasion_map.svg"):
    a_vals = np.linspace(0, 1, 100)

    # precompute resident SS for each species as function of its own a
    ss1 = []
    ss2 = []
    for a in a_vals:
        sp1.a = float(a)
        ss1.append(resident_ss(sp1, D, mode=UPTAKE_MODE))  # (s,g,alive,reason)

        sp2.a = float(a)
        ss2.append(resident_ss(sp2, D, mode=UPTAKE_MODE))

    Z = np.zeros((len(a_vals), len(a_vals)), dtype=int)  # rows=y=a_sp2, cols=x=a_sp1

    # codes:
    # 0 = neither invades (both alive)
    # 1 = sp2 invades sp1
    # 2 = sp1 invades sp2
    # 3 = mutual invasion (coexistence likely)
    # 4 = sp1 washes out (cannot be resident)
    # 5 = sp2 washes out
    # 6 = both wash out
    for iy, a2 in enumerate(a_vals):
        sp2.a = float(a2)
        s2, g2, alive2, _ = ss2[iy]

        for ix, a1 in enumerate(a_vals):
            sp1.a = float(a1)
            s1, g1, alive1, _ = ss1[ix]

            if (not alive1) and (not alive2):
                Z[iy, ix] = 6
                continue
            if not alive1:
                Z[iy, ix] = 4
                continue
            if not alive2:
                Z[iy, ix] = 5
                continue

            inv_2_into_1 = mu(sp2, s1, g1) > D
            inv_1_into_2 = mu(sp1, s2, g2) > D
            Z[iy, ix] = (1 if inv_2_into_1 else 0) + (2 if inv_1_into_2 else 0)

    labels = {
        0: "no invasion",
        1: f"{sp2.name} invades {sp1.name}",
        2: f"{sp1.name} invades {sp2.name}",
        3: "mutual invasion",
        4: f"{sp1.name} washout",
        5: f"{sp2.name} washout",
        6: "both washout",
    }

    colorscale = [
        [0.00, "#f2f2f2"],
        [0.14, "#f2f2f2"],  # 0
        [0.15, "#99bbff"],
        [0.29, "#99bbff"],  # 1
        [0.30, "#ff9999"],
        [0.44, "#ff9999"],  # 2
        [0.45, "#b6f2c2"],
        [0.59, "#b6f2c2"],  # 3
        [0.60, "#d9d2ff"],
        [0.74, "#d9d2ff"],  # 4
        [0.75, "#ffe0b3"],
        [0.89, "#ffe0b3"],  # 5
        [0.90, "#dddddd"],
        [1.00, "#dddddd"],  # 6
    ]

    fig = go.Figure(
        go.Heatmap(
            x=a_vals,
            y=a_vals,
            z=Z,
            zmin=0,
            zmax=6,
            colorscale=colorscale,
            colorbar=dict(
                tickmode="array",
                tickvals=list(labels.keys()),
                ticktext=[labels[k] for k in labels.keys()],
            ),
            hovertemplate=f"a_{sp1.name}=%{{x:.2f}}<br>a_{sp2.name}=%{{y:.2f}}<br>code=%{{z}}<extra></extra>",
        )
    )

    fig.update_layout(
        xaxis_title=f"a ({sp1.name})",
        yaxis_title=f"a ({sp2.name})",
        title=f"Invasion map (mode={UPTAKE_MODE}), D={D}",
        width=560,
        height=520,
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

at.a = 0.3
oa.a = 0.8

plot_zngi_and_consumption(at, oa, out="tmp.svg")
plot_invasion_map(at, oa, out="invasion_map.svg")
