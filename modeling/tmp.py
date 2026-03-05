import numpy as np
import plotly.graph_objects as go
from plotly.colors import qualitative

from style import *  # assumes style_plot(fig, ...) exists

# ----------------------------
# Parameters (species_resource)
# ----------------------------
V1 = 0.2
V2 = 0.2

K11 = 1.0  # sp1, res1
K12 = 0.5  # sp1, res2
K21 = 1.6  # sp2, res1
K22 = 0.1  # sp2, res2

# Stoichiometric requirements (resource per biomass formed)
Q11 = 0.8 / 10  # sp1, res1
Q12 = 1.0 / 10  # sp1, res2
Q21 = 1.0 / 10  # sp2, res1
Q22 = 0.7 / 10  # sp2, res2

mu_target = 0.1
S1, S2 = 6.0, 6.0  # supply point

R1s = np.linspace(0.0, 10.0, 200)
R2s = np.linspace(0.0, 10.0, 200)
R1g, R2g = np.meshgrid(R1s, R2s)

# ----------------------------
# Growth (SUBSTITUTABLE resources; additive Monod)
# IMPORTANT: PLUS (+) combines the two terms (no multiplication between resources)
# ----------------------------
mu1_r1 = V1 * (R1g / (K11 + R1g))
mu1_r2 = V1 * (R2g / (K12 + R2g))
mu1 = mu1_r1 + mu1_r2

mu2_r1 = V2 * (R1g / (K21 + R1g))
mu2_r2 = V2 * (R2g / (K22 + R2g))
mu2 = mu2_r1 + mu2_r2

# ----------------------------
# Consumption rates u_ij (species i, resource j)
# For substitutable resources, tie consumption of each resource to its own term.
# ----------------------------
u11 = Q11 * mu1_r1  # sp1 consumes res1
u12 = Q12 * mu1_r2  # sp1 consumes res2

u21 = Q21 * mu2_r1  # sp2 consumes res1
u22 = Q22 * mu2_r2  # sp2 consumes res2

# ----------------------------
# Consumption line through supply point
# If you want a Tilman-style constant "impact direction", use Q-ratios:
# slope = (ΔR2/ΔR1) = u12/u11 = Q12/Q11 (and similarly for species 2)
# ----------------------------
slope1 = Q12 / Q11
slope2 = Q22 / Q21

ys1 = S2 + slope1 * (R1s - S1)
ys2 = S2 + slope2 * (R1s - S1)


# Clip lines to plotting window (optional, but keeps it tidy)
def clip_line(y, ymin=0.0, ymax=10.0):
    y = y.copy()
    y[(y < ymin) | (y > ymax)] = np.nan
    return y


ys1 = clip_line(ys1, 0.0, R2s.max())
ys2 = clip_line(ys2, 0.0, R2s.max())

# ----------------------------
# Plot
# ----------------------------
fig = go.Figure()

# Consumption lines
for i, (y, slope) in enumerate([(ys1, slope1), (ys2, slope2)]):
    fig.add_trace(
        go.Scatter(
            x=R1s,
            y=y,
            mode="lines",
            line=dict(width=2, color=qualitative.Set1[i], dash="dash"),
            name=f"Species {i+1} consumption line (slope={slope:.2f})",
            hoverinfo="skip",
        )
    )

# ZNGIs: μ = mu_target
species = ["Species 1", "Species 2"]
for i, Z in enumerate([mu1, mu2]):
    fig.add_trace(
        go.Contour(
            x=R1s,
            y=R2s,
            z=Z,
            contours=dict(
                type="constraint",
                operation="=",
                value=mu_target,
                coloring="none",
                showlines=True,
            ),
            line=dict(width=2, color=qualitative.Set1[i]),
            showscale=False,
            name=f"{species[i]} ZNGI (μ={mu_target})",
            hoverinfo="skip",
        )
    )

# Supply point
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

fig.update_layout(
    title=f"Substitutable resources (additive Monod): ZNGIs μ={mu_target} 1/h",
    xaxis_title="Resource 1 (R1)",
    yaxis_title="Resource 2 (R2)",
    xaxis=dict(range=[0, R1s.max()]),
    yaxis=dict(range=[0, R2s.max()]),
    width=600,
    height=420,
)

# Tilman-style equal scaling
fig.update_yaxes(scaleanchor="x", scaleratio=1)

fig = style_plot(fig, line_thickness=3, marker_size=10)
fig.write_image("tmp.svg")
