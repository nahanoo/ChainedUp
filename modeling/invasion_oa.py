import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import plotly.graph_objects as go

from style import *


D = 0.3
PARAMETER_FILE = "parameters/parameters_15_mM_C_area.csv"


def load_params():
    return pd.read_csv(PARAMETER_FILE, index_col=0)


def make_resident_oa(params, a):
    p = params.loc["Oa"]
    return {
        "v_succ": float(p["v_succ"]),
        "v_gluc": float(p["v_gluc"]),
        "K_succ": float(p["K_succ"]),
        "K_gluc": float(p["K_gluc"]),
        "q_succ": float(p["q_succ"]),
        "q_gluc": float(p["q_gluc"]),
        "N0": float(p["N0"]),
        "s_in": float(p["s_in"]),
        "g_in": float(p["g_in"]),
        "a": float(a),
    }


def make_dummy_invader():
    return {
        "v_gluc": 0.0,
        "K_gluc": 1.0,
        "q_gluc": 1.0,
    }


def rhs(t, y, resident, invader):
    oa, inv, S, G = y

    a = resident["a"]

    J_oa_succ = a * resident["v_succ"] * S / (resident["K_succ"] + S)
    J_oa_gluc = (1.0 - a) * resident["v_gluc"] * G / (resident["K_gluc"] + G)
    mu_oa = J_oa_succ + J_oa_gluc

    mu_inv = invader["v_gluc"] * G / (invader["K_gluc"] + G)

    doa_dt = oa * (mu_oa - D)
    dinv_dt = inv * (mu_inv - D)

    dS_dt = D * (resident["s_in"] - S) - oa * (J_oa_succ / resident["q_succ"])
    dG_dt = (
        D * (resident["g_in"] - G)
        - oa * (J_oa_gluc / resident["q_gluc"])
        - inv * (mu_inv / invader["q_gluc"])
    )

    return [doa_dt, dinv_dt, dS_dt, dG_dt]


def simulate_resident_only(resident, t_end=200, n=2000):
    y0 = np.array(
        [
            resident["N0"],
            0.0,
            resident["s_in"],
            resident["g_in"],
        ],
        dtype=float,
    )

    t_eval = np.linspace(0, t_end, n)
    sol = solve_ivp(
        lambda t, y: rhs(t, y, resident, make_dummy_invader()),
        (0, t_end),
        y0,
        t_eval=t_eval,
        method="BDF",
        rtol=1e-8,
        atol=1e-12,
    )
    return sol.t, sol.y


def simulate_invasion(resident, invader, ti=50, inv_inoc=1e6, t_end=200, n=2000):
    t_eval1 = np.linspace(0, ti, max(2, int(n * ti / t_end)))
    t_eval2 = np.linspace(ti, t_end, max(2, int(n * (t_end - ti) / t_end)))

    y0_pre = np.array(
        [
            resident["N0"],
            0.0,
            resident["s_in"],
            resident["g_in"],
        ],
        dtype=float,
    )

    sol1 = solve_ivp(
        lambda t, y: rhs(t, y, resident, invader),
        (0, ti),
        y0_pre,
        t_eval=t_eval1,
        method="BDF",
        rtol=1e-8,
        atol=1e-12,
    )

    y_ti = sol1.y[:, -1].copy()
    y_ti[1] = inv_inoc

    sol2 = solve_ivp(
        lambda t, y: rhs(t, y, resident, invader),
        (ti, t_end),
        y_ti,
        t_eval=t_eval2,
        method="BDF",
        rtol=1e-8,
        atol=1e-12,
    )

    t = np.concatenate([sol1.t, sol2.t[1:]])
    y = np.concatenate([sol1.y, sol2.y[:, 1:]], axis=1)
    return t, y


def make_invader(params, G_star_a0, G_star_a05):
    """
    Synthetic glucose specialist invader in the same peak-area unit system.
    It is tuned so that it invades at a=0 but not at a=0.5.
    """
    if not (G_star_a0 > G_star_a05):
        raise ValueError(
            f"Need G*(a=0) > G*(a=0.5), got {G_star_a0:.4g} <= {G_star_a05:.4g}"
        )

    # Put break-even strictly between the two resident glucose equilibria
    G_crit = 0.5 * (G_star_a0 + G_star_a05)

    # plausible specialist affinity
    K_gluc_inv = 0.5  # mM

    # solve v_gluc so mu_inv(G_crit) = D
    v_gluc_inv = D * (K_gluc_inv + G_crit) / G_crit

    # use same biomass scale as Oa
    q_gluc_inv = float(params.loc["Oa", "q_gluc"])

    return {
        "v_gluc": v_gluc_inv,
        "K_gluc": K_gluc_inv,
        "q_gluc": q_gluc_inv,
        "G_crit": G_crit,
    }


def plot_species(fig, t, y, name):
    y = y.copy()
    y[y <= 0] = np.nan
    fig.add_trace(
        go.Scatter(
            x=t,
            y=y,
            mode="lines",
            name=name,
            line=dict(color=colors[name]),
        )
    )
    fig.update_layout(
        yaxis=dict(type="log", exponentformat="power", range=[3, 11]),
        legend=dict(
            x=0.05,
            y=0.05,
            xanchor="left",
            yanchor="bottom",
            bgcolor="rgba(255,255,255,0.5)",
        ),
    )
    return fig


def plot_resources(fig, t, y, name):
    fig.add_trace(
        go.Scatter(
            x=t,
            y=y,
            mode="lines",
            name=name,
            line=dict(color=colors[name]),
        )
    )
    fig.update_layout(
        legend=dict(
            x=0.96,
            y=0.5,
            xanchor="right",
            yanchor="middle",
            bgcolor="rgba(255,255,255,0.5)",
        ),
        yaxis=dict(range=[-3, 1.5], type="log", exponentformat="power"),
    )
    return fig


def update_layout(fig, title, yaxis_title):
    fig.update_layout(
        title=title,
        xaxis=dict(title="Time (h)", range=[0, 100], dtick=50),
        yaxis=dict(title=yaxis_title),
        width=150,
        height=150,
        showlegend=False,
    )
    fig = style_plot(fig, line_thickness=3)
    return fig


def resident_species_resources(params):
    resident_a0 = make_resident_oa(params, a=0.0)
    resident_a05 = make_resident_oa(params, a=0.5)

    t, y = simulate_resident_only(resident_a0)
    oa, inv, S, G = y
    fig = go.Figure()
    fig = plot_resources(fig, t, S, "Succinate")
    fig = plot_resources(fig, t, G, "Glucose")
    fig = update_layout(fig, "Oa, a = 0", "Area")
    fig.write_image("plots/invasion/resident_resources_oa_a0.svg")

    t, y = simulate_resident_only(resident_a05)
    oa, inv, S, G = y
    fig = go.Figure()
    fig = plot_resources(fig, t, S, "Succinate")
    fig = plot_resources(fig, t, G, "Glucose")
    fig = update_layout(fig, "Oa, a = 0.5", "Area")
    fig.write_image("plots/invasion/resident_resources_oa_a05.svg")

    return float(y[3, -1]), float(simulate_resident_only(resident_a0)[1][3, -1])


def invasion_dynamics(params):
    resident_a0 = make_resident_oa(params, a=0.0)
    resident_a05 = make_resident_oa(params, a=0.5)

    _, y0 = simulate_resident_only(resident_a0)
    _, y05 = simulate_resident_only(resident_a05)

    G_star_a0 = float(y0[3, -1])
    G_star_a05 = float(y05[3, -1])

    print(f"G*(a=0.0):  {G_star_a0:.6g}")
    print(f"G*(a=0.5):  {G_star_a05:.6g}")

    invader = make_invader(params, G_star_a0, G_star_a05)

    print("Synthetic invader parameters:")
    print(f"  v_gluc = {invader['v_gluc']:.6g}")
    print(f"  K_gluc = {invader['K_gluc']:.6g}")
    print(f"  q_gluc = {invader['q_gluc']:.6g}")
    print(f"  G_crit = {invader['G_crit']:.6g}")

    # invasion should succeed
    t, y = simulate_invasion(resident_a0, invader, ti=50, inv_inoc=1e6)
    oa, inv, S, G = y
    fig = go.Figure()
    fig = plot_species(fig, t, oa, "Oa")
    fig = plot_species(fig, t, inv, "At")
    fig = update_layout(fig, "Invasion dynamics", "CFUs/mL")
    fig.write_image("plots/invasion/invasion_dynamics_oa_a0.svg")

    # invasion should fail
    t, y = simulate_invasion(resident_a05, invader, ti=50, inv_inoc=1e6)
    oa, inv, S, G = y
    fig = go.Figure()
    fig = plot_species(fig, t, oa, "Oa")
    fig = plot_species(fig, t, inv, "At")
    fig = update_layout(fig, "Invasion dynamics", "CFUs/mL")
    fig.write_image("plots/invasion/invasion_dynamics_oa_a05.svg")


def main():
    params = load_params()
    invasion_dynamics(params)
    resident_species_resources(params)


if __name__ == "__main__":
    main()
