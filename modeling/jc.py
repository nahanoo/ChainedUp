import os
import numpy as np
import plotly.graph_objects as go
from plotly.colors import qualitative
from scipy.integrate import solve_ivp

from style import *  # style_plot(fig, ...)

# ----------------------------
# Parameters (species_resource)
# ----------------------------
V1 = 0.2
V2 = 0.2

K11 = 1.0  # sp1, res1
K12 = 0.5  # sp1, res2
K21 = 1.6  # sp2, res1
K22 = 0.1  # sp2, res2

# Re-parameterized stoichiometry: q_ij are coefficients (dimensionless)
q11 = 0.8  # sp1, res1
q12 = 1.0  # sp1, res2
q21 = 1.0  # sp2, res1
q22 = 0.7  # sp2, res2

# Yield: OD per mM
Y = 1.0  # OD / mM  (so 1/Y converts OD-growth to mM-consumption)

mu_target = 0.1  # use as chemostat loss rate D
D = mu_target

S1, S2 = 6.0, 6.0  # supply point / inflow concentrations (mM)

R1s = np.linspace(0.01, 10, 100)
R2s = np.linspace(0.01, 10, 100)
R1g, R2g = np.meshgrid(R1s, R2s)


# ----------------------------
# Growth (essential resources; multiplicative Monod)
# ----------------------------
def mu1_scalar(R1, R2):
    R1 = max(float(R1), 0.0)
    R2 = max(float(R2), 0.0)
    return V1 * (R1 / (K11 + R1)) * (R2 / (K12 + R2))


def mu2_scalar(R1, R2):
    R1 = max(float(R1), 0.0)
    R2 = max(float(R2), 0.0)
    return V2 * (R1 / (K21 + R1)) * (R2 / (K22 + R2))


# Grid versions for plotting ZNGIs
mu1 = V1 * (R1g / (K11 + R1g)) * (R2g / (K12 + R2g))
mu2 = V2 * (R1g / (K21 + R1g)) * (R2g / (K22 + R2g))

# ----------------------------
# Consumption vector lines through supply point
# (fixed impact direction due to stoichiometry coefficients)
# ----------------------------
slope1 = q12 / q11
ys1 = S2 + slope1 * (R1s - S1)

slope2 = q22 / q21
ys2 = S2 + slope2 * (R1s - S1)


# ----------------------------
# Chemostat ODE system (state = [R1, R2, N1, N2])
# ----------------------------
def model(t, y):
    R1, R2, N1, N2 = y

    m1 = mu1_scalar(R1, R2)
    m2 = mu2_scalar(R1, R2)

    # Biomass (OD)
    dN1dt = N1 * (m1 - D)
    dN2dt = N2 * (m2 - D)

    # Resources (mM)
    dR1dt = D * (S1 - R1) - (q11 * m1 * N1 + q21 * m2 * N2) / Y
    dR2dt = D * (S2 - R2) - (q12 * m1 * N1 + q22 * m2 * N2) / Y

    return [dR1dt, dR2dt, dN1dt, dN2dt]


# ----------------------------
# Invasion simulation:
# resident-only for t_res, then add invader at t_res
# ----------------------------
def integrate_invasion(
    t_res=200.0,
    t_end=2000.0,
    dt=0.2,
    n1_0=0.01,
    n2_inv=0.01,
    r1_0=None,
    r2_0=None,
):
    if r1_0 is None:
        r1_0 = S1
    if r2_0 is None:
        r2_0 = S2

    # Phase 1: species 1 only
    t_eval1 = np.arange(0.0, t_res + 1e-12, dt)
    y0_1 = [r1_0, r2_0, n1_0, 0.0]
    sol1 = solve_ivp(model, [0.0, t_res], y0_1, t_eval=t_eval1, rtol=1e-7, atol=1e-9)

    # State at invasion
    R1_res, R2_res, N1_res, _ = sol1.y[:, -1]

    # Phase 2: add invader (species 2)
    t_eval2 = np.arange(t_res, t_end + 1e-12, dt)
    y0_2 = [R1_res, R2_res, N1_res, n2_inv]
    sol2 = solve_ivp(model, [t_res, t_end], y0_2, t_eval=t_eval2, rtol=1e-7, atol=1e-9)

    # Stitch (avoid duplicate t_res point)
    if np.isclose(sol1.t[-1], sol2.t[0]):
        t_all = np.concatenate([sol1.t, sol2.t[1:]])
        y_all = np.concatenate([sol1.y, sol2.y[:, 1:]], axis=1)
    else:
        t_all = np.concatenate([sol1.t, sol2.t])
        y_all = np.concatenate([sol1.y, sol2.y], axis=1)

    return t_all, y_all  # y_all rows: R1, R2, N1, N2


def plot_invasion(out_dir="jc_plots", fname="invasion.pdf", t_res=200.0):
    os.makedirs(out_dir, exist_ok=True)
    t, y = integrate_invasion(t_res=t_res)
    R1, R2, N1, N2 = y

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=t,
            y=N1,
            mode="lines",
            line=dict(width=2, color=qualitative.Set1[0]),
            name="Species 1",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t,
            y=N2,
            mode="lines",
            line=dict(width=2, color=qualitative.Set1[1]),
            name="Species 2",
        )
    )

    # Mark invasion time
    fig.add_vline(
        x=t_res,
        line=dict(width=2, dash="dash", color="rgba(80,80,80,0.8)"),
        annotation_text="Invasion",
        annotation_position="top left",
    )

    fig.update_layout(
        title="Invasion in chemostat (Species 2 introduced into Species 1 steady state)",
        xaxis_title="Time (h)",
        yaxis_title="Biomass (OD)",
        width=600,
        height=400,
    )

    fig = style_plot(fig, line_thickness=3)
    fig.write_image(os.path.join(out_dir, fname))
    return fig


# ----------------------------
# Existing trajectory plot (kept)
# ----------------------------
def integrate_model():
    y0 = [S1, S2, 0.01, 0.01]  # start with small inoculum
    t_eval = np.linspace(0, 2000, 10000)
    sol = solve_ivp(model, [t_eval[0], t_eval[-1]], y0, t_eval=t_eval)
    return sol


def plot_trajectory(out_dir="jc_plots"):
    os.makedirs(out_dir, exist_ok=True)
    sol = integrate_model()
    N1, N2 = sol.y[2], sol.y[3]
    fig = go.Figure()
    for i, N in enumerate([N1, N2]):
        fig.add_trace(
            go.Scatter(
                x=sol.t,
                y=N,
                mode="lines",
                line=dict(width=2, color=qualitative.Set1[i]),
                name=f"Species {i+1}",
            )
        )
    fig.update_layout(
        title="Population Trajectories",
        xaxis_title="Time (h)",
        yaxis_title="Biomass (OD)",
        width=600,
        height=400,
    )
    fig = style_plot(fig, line_thickness=3)
    fig.write_image(os.path.join(out_dir, "trajectory.pdf"))
    return fig


# Run plots
plot_trajectory()
plot_invasion(fname="invasion.pdf", t_res=300.0)


# ----------------------------
# Plot helpers (ZNGIs) (unchanged)
# ----------------------------
def update_layout(fig):
    fig = style_plot(fig, line_thickness=3, marker_size=10)
    fig.update_layout(
        title=f"Zero Net Growth Isoclines (ZNGIs) μ = {mu_target} 1/h",
        xaxis_title="Resource 1 (R1)",
        yaxis_title="Resource 2 (R2)",
        xaxis=dict(range=[0, R1s.max()], dtick=2),
        yaxis=dict(range=[0, R2s.max()], dtick=2),
        width=600,
        height=400,
    )
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


def zngis(out_dir="jc_plots"):
    os.makedirs(out_dir, exist_ok=True)

    fig = go.Figure()

    fig.add_trace(
        go.Contour(
            x=R1s,
            y=R2s,
            z=mu1,
            contours=dict(
                type="constraint",
                operation="=",
                value=mu_target,
                coloring="none",
                showlines=True,
            ),
            line=dict(width=2, color=qualitative.Set1[0]),
            showscale=False,
            name="Species 1 ZNGI",
            hoverinfo="skip",
            showlegend=True,
        )
    )
    fig = update_layout(fig)
    fig.update_layout(width=500)
    fig.write_image(os.path.join(out_dir, "fig1.pdf"))

    fig.add_trace(
        go.Scatter(
            x=[S1],
            y=[S2],
            mode="markers+text",
            marker=dict(size=10, color=qualitative.Set1[2]),
            text=["Supply Point"],
            textposition="top left",
            showlegend=False,
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=R1s,
            y=ys1,
            mode="lines",
            line=dict(width=2, color=qualitative.Set1[0], dash="dash"),
            name="Species 1 Consumption Vector",
            hoverinfo="skip",
        )
    )
    fig = update_layout(fig)
    fig.write_image(os.path.join(out_dir, "fig2.pdf"))

    fig.add_trace(
        go.Contour(
            x=R1s,
            y=R2s,
            z=mu2,
            contours=dict(
                type="constraint",
                operation="=",
                value=mu_target,
                coloring="none",
                showlines=True,
            ),
            line=dict(width=2, color=qualitative.Set1[1]),
            showscale=False,
            name="Species 2 ZNGI",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=R1s,
            y=ys2,
            mode="lines",
            line=dict(width=2, color=qualitative.Set1[1], dash="dash"),
            name="Species 2 Consumption Vector",
            hoverinfo="skip",
        )
    )
    fig = update_layout(fig)
    fig.write_image(os.path.join(out_dir, "fig3.pdf"))

    return fig


# fig = zngis()
