from model_chain import Model
from experiment import Species
import numpy as np
from style import *
import plotly.graph_objects as go
import pandas as pd
from os import path
from multiprocessing import Pool
from plotly.subplots import make_subplots
from scipy.optimize import least_squares

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

    feed_g = float(params["g_in"].iloc[0])
    feed_s = float(params["s_in"].iloc[0])
    model.M_glucose = feed_g
    model.M_succinate = feed_s

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


def pick_evenly_spaced_indices(x, y, n_points=50, log_space=True, eps=1e-12):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) == 0:
        return np.array([], dtype=int)

    if log_space:
        X = np.log10(np.clip(x, eps, None))
        Y = np.log10(np.clip(y, eps, None))
    else:
        X = x.copy()
        Y = y.copy()

    keep = np.ones(len(X), dtype=bool)
    keep[1:] = (np.diff(X) != 0) | (np.diff(Y) != 0)
    X = X[keep]
    Y = Y[keep]
    original_idx = np.where(mask)[0][keep]

    if len(original_idx) <= n_points:
        return original_idx

    ds = np.sqrt(np.diff(X) ** 2 + np.diff(Y) ** 2)
    s = np.concatenate([[0.0], np.cumsum(ds)])

    if s[-1] == 0:
        idx = np.linspace(0, len(original_idx) - 1, n_points).astype(int)
        return original_idx[idx]

    s_targets = np.linspace(0.0, s[-1], n_points)
    idx = np.searchsorted(s, s_targets, side="left")
    idx = np.clip(idx, 0, len(original_idx) - 1)
    idx = np.unique(idx)

    return original_idx[idx]


def a_mono_culture():
    conc = 15
    p_f = path.join("parameters", "parameters_15_mM_C_area.csv")
    params = pd.read_csv(p_f, index_col=0)

    xs = np.linspace(0, 100, 1000)

    # simulate densely in a
    aas_fine = np.linspace(0.0, 1.0, 2000)
    i_half = np.argmin(np.abs(aas_fine - 0.5))
    with Pool() as pool:
        args_oa = [(conc, params, xs, "Oa", a) for a in aas_fine]
        results_oa = np.asarray(pool.map(simulate_endpoints, args_oa))

    suc_oa_c1, gluc_oa_c1, _, oa_c1, suc_oa_c2, gluc_oa_c2, _, oa_c2 = results_oa.T

    # keep only evenly spaced points in displayed log-log space
    idx_oa = pick_evenly_spaced_indices(
        suc_oa_c1,
        gluc_oa_c1,
        n_points=30,
        log_space=True,
    )

    eps = 1e-12

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=np.clip(suc_oa_c1[idx_oa], eps, None),
            y=np.clip(gluc_oa_c1[idx_oa], eps, None),
            mode="markers",
            marker=dict(
                color=aas_fine[idx_oa],
                colorscale="cividis",
                showscale=True,
                size=9,
                colorbar=dict(
                    title="a",
                    len=0.8,
                    lenmode="fraction",
                    y=0.5,
                    yanchor="middle",
                    thickness=14,
                ),
            ),
            name="Oa",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[np.clip(suc_oa_c1[i_half], eps, None)],
            y=[np.clip(gluc_oa_c1[i_half], eps, None)],
            mode="markers",
            marker=dict(
                size=16,
                symbol="circle-open",
                line=dict(width=3, color="red"),
            ),
            name="a = 0.5",
            showlegend=True,
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Succinate (Area)", type="log", exponentformat="power", range=[-3, 1]
        ),
        yaxis=dict(
            title="Glucose (Area)", type="log", exponentformat="power", range=[-3, 1]
        ),
        title="Steady states at different allocations",
        height=200,
        width=240,
        showlegend=False,
    )

    fig = style_plot(fig, line_thickness=3, marker_size=8, font_size=11)
    fig.write_image("plots/contours/mono_culture_resource_allocation.svg")


a_mono_culture()


def resource_allocation_heatmap():
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
    fig_mono = go.Figure()
    fig_mono.add_trace(
        go.Heatmap(
            x=aas,
            y=aas,
            z=z_c1,
            colorscale="cividis",
            colorbar=dict(
                title="log<sub>10</sub>(At/Oa)",
            ),
            showscale=False,
        ),
    )
    fig_mono.update_layout(
        xaxis=dict(title="Resource allocation At", dtick=0.2),
        xaxis2=dict(title="Resource allocation At", dtick=0.2),
        yaxis=dict(title="Resource allocation Oa", dtick=0.2),
        width=380,
        height=380,
        title="Coexistence in succinate-glucose space",
    )
    fig_mono = style_plot(fig_mono, font_size=14, line_thickness=2, top_margin=35)
    fig_mono.write_image("plots/contours/coexistence_resource_allocation_mono.svg")


def at_resource_allocation_steady_state():
    aas = np.linspace(0.0, 1.0, 20)
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    xs = np.linspace(0, 100, 1000)
    with Pool() as pool:
        args = [(conc, params, xs, "At", a_at) for a_at in aas]
        results = pool.map(simulate_endpoints, args)
        suc_c1, gluc_c1, at_c1, oa_c1, suc_c2, gluc_c2, at_c2, oa_c2 = np.asarray(
            results
        ).T

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=aas,
            y=at_c1,
            mode="lines+markers",
            marker=dict(color=colors["Succinate+Glucose"]),
            name="Chemostat 1",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=aas,
            y=at_c2,
            mode="lines+markers",
            marker=dict(color=colors["Succinate+Glucose Outflow"]),
            name="Downstream<br>Chemostat",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Resource allocation At"),
        yaxis=dict(title="At steady state"),
        title="At steady state as function of resource allocation",
        legend=dict(y=0.55, x=0.5, yanchor="top", xanchor="center"),
        height=190,
        width=380,
    )
    fig = style_plot(fig, font_size=12, marker_size=8, line_thickness=2)
    fig.write_image("plots/contours/at_resource_allocation.svg")


def carbon_allocation():
    aas = np.linspace(0.0, 1.0, 20)
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    xs = np.linspace(0, 100, 1000)
    with Pool() as pool:
        args = [(conc, params, xs, "At", a_at) for a_at in aas]
        results = pool.map(simulate_endpoints, args)
        suc_c1, gluc_c1, at_c1, oa_c1, suc_c2, gluc_c2, at_c2, oa_c2 = np.asarray(
            results
        ).T

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=aas,
            y=4 * suc_c1,
            mode="lines+markers",
            marker=dict(color=colors["Succinate"]),
            name="Carbon from Succinate",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=aas,
            y=6 * gluc_c1,
            mode="lines+markers",
            marker=dict(color=colors["Glucose"]),
            name="Carbon from Glucose",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=aas,
            y=6 * gluc_c1 + 4 * suc_c1,
            mode="lines+markers",
            marker=dict(color=colors["Succinate+Glucose"]),
            name="Total carbon",
        )
    )
    fig.update_layout(
        xaxis=dict(title="Resource allocation At"),
        yaxis=dict(title="Total carbon in chemostat 1 [mM C]"),
        title="Total carbon in chemostat 1 as function of resource allocation",
    )
    fig = style_plot(fig, font_size=12, marker_size=8, line_thickness=2)
    fig.write_image("plots/contours/at_total_carbon.svg")


def effective_allocation():
    aas = np.linspace(0.0, 1.0, 20)
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    xs = np.linspace(0, 100, 1000)
    with Pool() as pool:
        args = [(conc, params, xs, "At", a_at) for a_at in aas]
        results = pool.map(simulate_endpoints, args)
        suc_c1, gluc_c1, at_c1, oa_c1, suc_c2, gluc_c2, at_c2, oa_c2 = np.asarray(
            results
        ).T
    atp = params.loc["At"]
    v_succ = float(atp["v_succ"])
    K_succ = float(atp["K_succ"])
    v_gluc = float(atp["v_gluc"])
    K_gluc = float(atp["K_gluc"])

    # fluxes at endpoints
    J_S_1 = aas * v_succ * suc_c1 / (K_succ + suc_c1)
    J_G_1 = (1 - aas) * v_gluc * gluc_c1 / (K_gluc + gluc_c1)

    J_S_2 = aas * v_succ * suc_c2 / (K_succ + suc_c2)
    J_G_2 = (1 - aas) * v_gluc * gluc_c2 / (K_gluc + gluc_c2)

    a_eff_1 = J_S_1 / (J_S_1 + J_G_1)
    a_eff_2 = J_S_2 / (J_S_2 + J_G_2)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=aas,
            y=a_eff_1,
            mode="lines+markers",
            name="a_eff in C1",
            marker=dict(color=colors["Succinate"]),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=aas,
            y=a_eff_2,
            mode="lines+markers",
            name="a_eff in C2",
            marker=dict(color=colors["Glucose"]),
        )
    )
    fig.update_layout(
        xaxis_title="parameter a (At)",
        yaxis_title="effective allocation a_eff = J_S / (J_S + J_G)",
        title="Allocation parameter vs realized uptake split",
    )
    fig = style_plot(fig, font_size=12, marker_size=8, line_thickness=2)
    fig.write_image("plots/contours/at_effective_allocation.svg")


def compute_fluxes(S_grid, G_grid, sp):
    return sp.a * sp.v_succ * S_grid / (sp.K_succ + S_grid) + (
        1 - sp.a
    ) * sp.v_gluc * G_grid / (sp.K_gluc + G_grid)


def feed_from_conc(conc):
    # same mapping as your Model class
    C_to_mM_succinate = {45: 11.25, 30: 7.5, 15: 3.75, 7.5: 1.875, 5: 1.25, 2.5: 0.625}
    C_to_mM_glucose = {45: 7.5, 30: 5, 15: 2.5, 7.5: 1.25, 5: 0.833, 2.5: 0.417}
    return C_to_mM_succinate[conc], C_to_mM_glucose[conc]  # (S0, G0)


def steady_state_Rstar_fixed_a(sp, *, a, D, S0, G0):
    """
    Compute (S*, G*) for a monoculture chemostat at fixed allocation a.

    Model:
      J_S = a*v_succ*s/(K_succ+s)
      J_G = (1-a)*v_gluc*g/(K_gluc+g)
      J = J_S + J_G

    Steady state with N>0:
      (1) J(S*,G*) = D
      (2) N inferred from glucose balance equals N inferred from succinate balance:
          D(G0-G*) q_gluc / J_G  =  D(S0-S*) q_succ / J_S
    """
    sp.a = a

    J_max = sp.a * sp.v_succ + (1.0 - sp.a) * sp.v_gluc
    if D >= J_max:
        raise ValueError(f"Washout: D={D} >= J_max={J_max:.3g} for a={a}")

    def J_S(s):
        return sp.a * sp.v_succ * s / (sp.K_succ + s)

    def J_G(g):
        return (1.0 - sp.a) * sp.v_gluc * g / (sp.K_gluc + g)

    def residual(x):
        s, g = x  # NOTE: order matches plot axes (x=succinate, y=glucose)
        js = J_S(s)
        jg = J_G(g)
        j = js + jg

        r1 = j - D

        eps = 1e-12
        # N* from each resource balance (same as in your chemostat equations)
        N_from_g = D * (G0 - g) * sp.q_gluc / (jg + eps)
        N_from_s = D * (S0 - s) * sp.q_succ / (js + eps)
        r2 = N_from_g - N_from_s

        return np.array([r1, r2], dtype=float)

    # initial guess: small succinate (since Km is tiny) and moderate glucose
    x0 = np.array([min(S0, 1e-2), 0.5 * G0], dtype=float)

    res = least_squares(
        residual,
        x0=x0,
        bounds=([0.0, 0.0], [S0, G0]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    s_star, g_star = res.x
    return float(s_star), float(g_star)


def plot_isoclines():
    conc = 15
    level = 0.3  # this is your ZNGI/isocline level; interpret as D (1/h) for R*
    D = level
    allocation = 0.5

    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    oa = Species("Oa", params.loc["Oa"])
    sp = oa

    S0, G0 = 10, 1.5

    names = ["Glucose", "Succinate+Glucose", "Succinate"]
    a_values = [0.0, allocation, 1.0]

    fig = go.Figure()

    # log grid for succinate to resolve Km ~ 1e-3
    S_concs = np.logspace(-6, np.log10(S0), 250)
    # keep glucose linear (you can switch to log similarly if you want)
    G_concs = np.linspace(0.0, G0, 250)

    S_grid, G_grid = np.meshgrid(S_concs, G_concs)

    for i, a in enumerate(a_values):
        name = names[i]
        sp.a = a
        flux_grid = compute_fluxes(S_grid, G_grid, sp)

        fig.add_trace(
            go.Contour(
                x=S_concs,
                y=G_concs,
                z=flux_grid,
                contours=dict(
                    type="constraint",
                    operation="=",
                    value=level,
                    coloring="none",
                    showlines=True,
                ),
                line=dict(width=2, color=colors[name]),
                showscale=False,
                name=f"{name}: a={a:.1f}",
                hoverinfo="skip",
            )
        )

    # ---- Add the chemostat steady-state point (R*) for a=0.5 only ----
    a_star = allocation
    s_star, g_star = steady_state_Rstar_fixed_a(sp, a=a_star, D=D, S0=S0, G0=G0)

    fig.add_trace(
        go.Scatter(
            x=[s_star],
            y=[g_star],
            mode="markers+text",
            marker=dict(size=11, color=colors["Succinate+Glucose"], line=dict(width=1)),
            text=[f"$S^*={s_star:.2g}, G^*={g_star:.2g}$"],
            textposition="top left",
            showlegend=False,
            hoverinfo="skip",
        )
    )

    # Optional: show the feed point for context
    fig.add_trace(
        go.Scatter(
            x=[S0],
            y=[G0],
            mode="markers+text",
            marker=dict(size=9, color="rgba(80,80,80,0.9)"),
            text=["Feed"],
            textposition="bottom right",
            showlegend=False,
            hoverinfo="skip",
        )
    )

    fig.update_layout(
        xaxis=dict(title="Succinate [mM]", type="log", exponentformat="power"),
        yaxis=dict(title="Glucose [mM]"),
        legend=dict(
            yref="paper",
            xref="paper",
            x=0.99,
            y=0.99,
            xanchor="right",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.5)",
        ),
        width=400,
        height=360,
    )

    fig = style_plot(fig, font_size=font_size, line_thickness=3, marker_size=10)
    fig.write_image("plots/contours/isoclines_f_a.svg")
    return fig


def G_required_for_level(S, sp, level):
    """
    For each S, solve for G such that:
        a*v_succ*S/(K_succ+S) + (1-a)*v_gluc*G/(K_gluc+G) = level

    Returns G(S) (NaN where infeasible).
    """
    S = np.asarray(S, dtype=float)
    G = np.full_like(S, np.nan, dtype=float)

    # succinate contribution
    J_S = sp.a * sp.v_succ * S / (sp.K_succ + S)

    # remaining needed from glucose
    rem = level - J_S

    # If rem <= 0, succinate alone already achieves the level (any G>=0 works)
    # Here we set G=0 to visualize "glucose can be driven to ~0".
    mask_zero = rem <= 0
    G[mask_zero] = 0.0

    # If rem >= max glucose contribution, infeasible (no solution)
    JG_max = (1 - sp.a) * sp.v_gluc
    mask_bad = rem > JG_max
    # keep NaN

    # Solve Monod explicitly for the feasible region: rem in (0, JG_max]
    mask = (rem > 0) & (rem <= JG_max)

    # rem = (1-a)*v_gluc * G/(K+G)  ->  G = rem*K / (JG_max - rem)
    K = sp.K_gluc
    G[mask] = rem[mask] * K / (JG_max - rem[mask])

    return G


def plot_required_glucose_vs_succinate():
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)

    oa = Species("Oa", params.loc["Oa"])
    sp = oa

    level = 0.3
    a_values = [0.0, 0.5, 1.0]
    names = ["Glucose", "Succinate+Glucose", "Succinate"]

    # Use a log grid for succinate to resolve Km ~ 1e-3
    S = np.logspace(-6, np.log10(3.75), 400)

    fig = go.Figure()

    for i, a in enumerate(a_values):
        sp.a = a
        name = names[i]

        G_req = G_required_for_level(S, sp, level)

        fig.add_trace(
            go.Scatter(
                x=S,
                y=G_req,
                mode="lines",
                line=dict(width=2, color=colors[name]),
                name=f"{name}: a={a:.1f}",
                hoverinfo="skip",
            )
        )

    # Reference lines at Km
    fig.add_vline(
        x=sp.K_succ,
        line=dict(width=2, dash="dash", color="rgba(80,80,80,0.8)"),
        annotation_text="K_succ",
        annotation_position="top left",
    )
    fig.add_hline(
        y=sp.K_gluc,
        line=dict(width=2, dash="dash", color="rgba(80,80,80,0.8)"),
        annotation_text="K_gluc",
        annotation_position="bottom right",
    )

    fig.update_layout(
        title=f"Required glucose along isocline (J={level}) vs succinate",
        xaxis=dict(title="Succinate [mM]", type="log", exponentformat="power"),
        yaxis=dict(title="Glucose required [mM]"),
        legend=dict(
            yref="paper",
            xref="paper",
            x=0.99,
            y=0.99,
            xanchor="right",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.5)",
        ),
        width=200 * 2,
        height=180 * 2,
    )

    fig = style_plot(fig, font_size=font_size, line_thickness=3, marker_size=10)
    fig.write_image("plots/contours/required_glucose_vs_succinate.svg")
