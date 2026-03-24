import pandas as pd
import numpy as np
from style import *
import numpy as np
from scipy.integrate import solve_ivp
import plotly.graph_objects as go

V11 = 0.6
V12 = 0.5
K11 = 0.3
K12 = 0.5
Q11 = 0.1
Q12 = 0.1

# invader (species 2) - glucose specialist
V21 = 0.65
K21 = 0.5
Q21 = 0.1

# chemostat
D = 0.3
S0 = 5.0  # give plenty of the "other" resource
G0 = 2  # MUST be > 1.5 for the glucose specialist to survive


def rhs(t, y, a):
    N1, N2, S, G = y

    J11 = a * V11 * S / (K11 + S)
    J12 = (1 - a) * V12 * G / (K12 + G)
    mu1 = J11 + J12

    mu2 = V21 * G / (K21 + G)

    dN1dt = N1 * (mu1 - D)
    dN2dt = N2 * (mu2 - D)

    dSdt = D * (S0 - S) - N1 * (J11 / Q11)
    dGdt = D * (G0 - G) - N1 * (J12 / Q12) - N2 * (mu2 / Q21)

    return [dN1dt, dN2dt, dSdt, dGdt]


def simulate(t0, ti, y0, a, N2_inoc=0.1, n=2000, t_end=200):
    t_eval1 = np.linspace(t0, ti, max(2, int(n * (ti - t0) / (t_end - t0))))
    t_eval2 = np.linspace(ti, t_end, max(2, int(n * (t_end - ti) / (t_end - t0))))

    y0_pre = y0.copy()
    y0_pre[1] = 0.0  # N2 = 0
    sol1 = solve_ivp(lambda t, y: rhs(t, y, a), (t0, ti), y0_pre, t_eval=t_eval1)

    y_ti = sol1.y[:, -1].copy()
    y_ti[1] = N2_inoc

    sol2 = solve_ivp(lambda t, y: rhs(t, y, a), (ti, t_end), y_ti, t_eval=t_eval2)

    t = np.concatenate([sol1.t, sol2.t[1:]])
    y = np.concatenate([sol1.y, sol2.y[:, 1:]], axis=1)
    return t, y


def plot_species(fig, t, y, name):
    y[y == 0] = np.nan
    fig.add_trace(
        go.Scatter(x=t, y=y, mode="lines", name=name, line=dict(color=colors[name]))
    )
    fig.update_layout(
        yaxis=dict(
            range=[-3, 0],
            type="log",
        ),
        legend=dict(
            x=0.05,
            y=0.05,
            xanchor="left",
            yanchor="bottom",
            bgcolor="rgba(255, 255, 255, 0.5)",
        ),
    )
    return fig


def plot_resources(fig, t, y, name):
    fig.add_trace(
        go.Scatter(x=t, y=y, mode="lines", name=name, line=dict(color=colors[name]))
    )
    fig.update_layout(
        yaxis=dict(dtick=1),
        legend=dict(
            x=0.96,
            y=0.5,
            xanchor="right",
            yanchor="middle",
            bgcolor="rgba(255, 255, 255, 0.5)",
        ),
    )
    return fig


def update_layout(fig, title, yaxis_title):
    fig.update_layout(
        title=title,
        xaxis=dict(title="Time (h)", range=[0, 100], dtick=50),
        yaxis=dict(title=yaxis_title),
        width=145,
        height=135,
    )
    fig = style_plot(fig, line_thickness=3)
    return fig


def resident_species_resources():
    y0 = np.array([0.05, 0.0, S0, G0])
    t, y = simulate(t0=0, ti=10, t_end=200, y0=y0, a=0.5, N2_inoc=0)
    N1, N2, S, G = y
    fig = go.Figure()
    fig = plot_resources(fig, t, S, "Succinate")
    fig = plot_resources(fig, t, G, "Glucose")
    fig = update_layout(fig, "Resident resources (generalist)", "Concentration (mM)")
    fig.write_image("plots/invasion/resident_resources_generalist.svg")

    t, y = simulate(t0=0, ti=10, t_end=200, y0=y0, a=0, N2_inoc=0)
    N1, N2, S, G = y
    fig = go.Figure()
    fig = plot_resources(fig, t, S, "Succinate")
    fig = plot_resources(fig, t, G, "Glucose")
    fig = update_layout(fig, "Resident resources (specialist)", "Concentration (mM)")
    fig.write_image("plots/invasion/resident_resources_specialist.svg")

    t, y = simulate(t0=0, ti=10, t_end=200, y0=y0, a=1, N2_inoc=0)
    N1, N2, S, G = y
    fig = go.Figure()
    fig = plot_resources(fig, t, S, "Succinate")
    fig = plot_resources(fig, t, G, "Glucose")
    fig = update_layout(fig, "Resident resources (specialist)", "Concentration (mM)")
    fig.write_image("plots/invasion/resident_resources_specialist_succinate.svg")


def invasion_dynamics():
    y0 = np.array([0.05, 0.0, S0, G0])
    t, y = simulate(t0=0, ti=50, t_end=200, y0=y0, a=0.5, N2_inoc=0.01)
    N1, N2, S, G = y
    fig = go.Figure()
    fig = plot_species(fig, t, N1, "Oa")
    fig = plot_species(fig, t, N2, "At")
    fig = update_layout(
        fig, "Invasion dynamics (generalist resident)", "Abundance (OD)"
    )
    fig.write_image("plots/invasion/invasion_dynamics_generalist.svg")

    t, y = simulate(t0=0, ti=50, t_end=200, y0=y0, a=0, N2_inoc=0.01)
    N1, N2, S, G = y
    fig = go.Figure()
    fig = plot_species(fig, t, N1, "Oa")
    fig = plot_species(fig, t, N2, "At")
    fig = update_layout(
        fig, "Invasion dynamics (specialist resident)", "Abundance (OD)"
    )
    fig.write_image("plots/invasion/invasion_dynamics_specialist.svg")


resident_species_resources()
invasion_dynamics()
