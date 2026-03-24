import numpy as np
import plotly.graph_objects as go
import pandas as pd
from os import path
from scipy.integrate import solve_ivp

from experiment import Species
from style import *


# ----------------------------
# Load parameters
# ----------------------------
conc = 15
p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)

D = 0.3

# Feed concentrations (mM) for this conc (your mapping)
C_to_mM_succinate = {45: 11.25, 30: 7.5, 15: 3.75, 7.5: 1.875, 5: 1.25, 2.5: 0.625}
C_to_mM_glucose = {45: 7.5, 30: 5, 15: 2.5, 7.5: 1.25, 5: 0.833, 2.5: 0.417}

G_in = C_to_mM_glucose[conc]
S_in = C_to_mM_succinate[conc]


# ----------------------------
# Monoculture Model (At only) with switch in allocation a(t)
# State: y = [N, g, s]
# ----------------------------
def simulate_at_switch(at: Species, t_switch: float, a_pre: float, a_post: float):
    def rhs(t, y):
        N, g, s = y

        a = a_pre if t < t_switch else a_post

        J_S = a * at.v_succ * s / (at.K_succ + s)
        J_G = (1.0 - a) * at.v_gluc * g / (at.K_gluc + g)
        J = J_S + J_G

        dNdt = N * (J - D)

        dgdt = D * (G_in - g) - (J_G * N / at.q_gluc)
        dsdt = D * (S_in - s) - (J_S * N / at.q_succ)

        return [dNdt, dgdt, dsdt]

    return rhs


# ----------------------------
# Run simulation
# ----------------------------
at = Species("At", params.loc["At"])
oa = Species(
    "Oa", params.loc["Oa"]
)  # not used in this simulation, but could be for later

t_switch = 50.0
t_end = 150.0  # change as needed
t = np.linspace(0, t_end, 1200)

# Initial conditions: start at feed concentrations, OD inoculum
N0 = 0.1
y0 = [N0, G_in, S_in]

rhs = simulate_at_switch(oa, t_switch=t_switch, a_pre=0.0, a_post=0.4)

sol = solve_ivp(
    rhs,
    t_span=(t[0], t[-1]),
    y0=y0,
    t_eval=t,
    method="BDF",
    rtol=1e-8,
    atol=1e-12,
)

N = sol.y[0]
g = sol.y[1]
s = sol.y[2]


# ----------------------------
# Plot (your style)
# ----------------------------
fig = go.Figure()

fig.add_trace(
    go.Scatter(
        x=sol.t,
        y=g,
        mode="lines",
        name="Glucose",
        line=dict(color=colors["Glucose"]),
    )
)
fig.add_trace(
    go.Scatter(
        x=sol.t,
        y=s,
        mode="lines",
        name="Succinate",
        line=dict(color=colors["Succinate"]),
    )
)
fig.add_vline(
    x=t_switch,
    line=dict(width=2, dash="dash", color="rgba(80,80,80,0.9)"),
    annotation_text="switch a: 0 → 0.5",
    annotation_position="top left",
)

fig.update_layout(
    title="At monoculture: allocation switch (glucose-only → dual substrate)",
    xaxis=dict(title="Time (h)"),
    yaxis=dict(title="Concentration (mM)"),
    legend=dict(x=0.96, y=0.5, xanchor="right", yanchor="middle"),
    width=450,
    height=320,
)

fig = style_plot(fig, line_thickness=2)
fig.write_image("plots/invasion/at_allocation_switch_resources.svg")


# Optional: plot biomass too
fig2 = go.Figure()
fig2.add_trace(
    go.Scatter(
        x=sol.t,
        y=N,
        mode="lines",
        name="At (OD)",
        line=dict(color=colors.get("At", qualitative.Set1[0])),
    )
)
fig2.add_vline(
    x=t_switch,
    line=dict(width=2, dash="dash", color="rgba(80,80,80,0.9)"),
    annotation_text="switch",
    annotation_position="top left",
)
fig2.update_layout(
    title="At monoculture biomass (OD)",
    xaxis=dict(title="Time (h)"),
    yaxis=dict(title="Biomass (OD)"),
    width=450,
    height=320,
)
fig2 = style_plot(fig2, line_thickness=2)
fig2.write_image("plots/invasion/at_allocation_switch_biomass.svg")
