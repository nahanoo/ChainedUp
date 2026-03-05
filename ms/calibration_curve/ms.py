from style import *
import plotly.graph_objects as go
import numpy as np

steps = range(16)
start = 2.5
current_element = start
y = []

for i in steps:
    y.append(current_element)
    current_element = current_element / 2

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=np.array(steps) + 1,
        y=y,
        mode="lines+markers",
        line=dict(color=colors["Glucose"], width=2),
        name="Glucose Concentration",
        showlegend=False,
    )
)
fig.add_trace(
    go.Scatter(
        x=np.array(list(steps)[::3]) + 1,
        y=y[::3],
        mode="markers",
        marker=dict(size=8, color=colors["Succinate"]),
        showlegend=False,
    )
)
fig.update_layout(
    xaxis_title="Step",
    yaxis_title="Glucose Concentration (mM)",
    yaxis=dict(type="log", tickformat=".1g", dtick=1),
)
fig = style_plot(fig, marker_size=6, line_thickness=2, font_size=14)
fig.write_image("plots/ms/gradient.svg")
