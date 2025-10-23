from page import model, parse_params
from scipy.integrate import odeint
import numpy as np
from style import *
import plotly.graph_objects as go

font_size = 14


def a_mono_culture():
    p = parse_params("parameters_dash_init.csv")
    p["D"] = 0.3
    t = np.linspace(0, 200, 2000)
    y0 = [0.1, 0.0, 1, 1, 0, 0, 0, 0]
    N1, N2, C1, C2, U1, U2, S1, S2 = odeint(model, y0, t, args=(p,)).T

    aas = np.linspace(0, 1, 100)
    C1s, C2s = [], []
    for a in aas:
        p["a1"] = a
        N1, N2, C1, C2, U1, U2, S1, S2 = odeint(model, y0, t, args=(p,)).T
        C1s.append(C1[-1])
        C2s.append(C2[-1])

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=C1s,
            y=C2s,
            mode="markers",
            marker=dict(
                symbol="triangle-up",
                color=aas,
                colorscale="viridis",
                colorbar=dict(
                    title="a",
                    len=0.55,  # height as fraction of plotting area
                    lenmode="fraction",
                    y=0.5,
                    yanchor="middle",
                    thickness=14,  # width in px (optional)
                ),
            ),
            name="At",
        )
    )
    y0 = [0.0, 0.1, 1, 1, 0, 0, 0, 0]
    N1, N2, C1, C2, U1, U2, S1, S2 = odeint(model, y0, t, args=(p,)).T

    C1s, C2s = [], []
    for a in aas:
        p["a2"] = a
        N1, N2, C1, C2, U1, U2, S1, S2 = odeint(model, y0, t, args=(p,)).T
        C1s.append(C1[-1])
        C2s.append(C2[-1])

    fig.add_trace(
        go.Scatter(
            x=C1s,
            y=C2s,
            mode="markers",
            marker=dict(
                color=aas,
                colorscale="viridis",
                showscale=False,
            ),
            line=dict(width=2),
            name="Oa",
        )
    )
    fig.update_layout(
        xaxis=dict(
            title="Succinate [mM]",
            # type="log"
        ),
        yaxis=dict(
            title="Glucose [mM]",
            # type="log"
        ),
        title="Mono-culture",
    )
    fig = style_plot(fig, line_thickness=3, marker_size=10, font_size=font_size)
    fig.write_image("plots/contours/a_C1_C2_resource_space.pdf")


def a1_a2_sweep():
    p = parse_params("parameters_dash_init.csv")
    p["D"] = 0.3
    y0 = [0.1, 0.1, 1, 1, 0, 0, 0, 0]
    t = np.linspace(0, 2000, 20000)
    a1_vals = np.linspace(0.0, 1.0, 100)
    a2_vals = np.linspace(0.0, 1.0, 100)

    def simulate_tail_means(a1, a2, p_base):
        p = dict(p_base)
        p["a1"] = float(a1)
        p["a2"] = float(a2)
        N1, N2, C1, C2, U1, U2, S1, S2 = odeint(model, y0, t, args=(p,), mxstep=5000).T
        i0 = int(0.95 * len(t))  # last 5%
        return (N1[i0:].mean(), N2[i0:].mean(), C1[i0:].mean(), C2[i0:].mean())

    # Sweep
    Z_C1 = np.empty((len(a2_vals), len(a1_vals)))
    Z_C2 = np.empty((len(a2_vals), len(a1_vals)))
    Z_N1 = np.empty_like(Z_C1)
    Z_N2 = np.empty_like(Z_C1)

    for i, a2 in enumerate(a2_vals):
        for j, a1 in enumerate(a1_vals):
            N1m, N2m, C1m, C2m = simulate_tail_means(a1, a2, p)
            Z_N1[i, j], Z_N2[i, j] = N1m, N2m
            Z_C1[i, j], Z_C2[i, j] = C1m, C2m

    # Coexistence mask (tolerance relative to initial biomass)
    tol = 1e-4  # adjust to your scale; e.g., 1% of initial -> 0.001 if y0 ~ 0.1
    coexist = (Z_N1 > tol) & (Z_N2 > tol)
    Zs = [Z_C1, Z_C2, Z_C1 / (Z_C1 + Z_C2)]
    scale_names = ["Succinate<br>[mM]", "Glucose<br>[mM]", "Succinate<br>ratio"]
    filenames = ["a1_a2_C1.pdf", "a1_a2_C2.pdf", "a1_a2_C1_C2_ratio.pdf"]
    # Your ratio plot
    for i, Z in enumerate(Zs):
        fig = go.Figure()
        fig.add_trace(
            go.Contour(
                z=Z,
                x=a1_vals,
                y=a2_vals,
                contours=dict(showlines=False),
                colorscale="Viridis",
                zmin=0,
                zmax=np.max(Zs[i]),
                ncontours=50,
                zauto=False,
                colorbar=dict(title=scale_names[i], len=0.6, y=0.2, thickness=12),
                name=scale_names[i],
            )
        )

        # Overlay isolines where each species crosses the persistence threshold
        fig.add_trace(
            go.Contour(
                z=Z_N1,
                x=a1_vals,
                y=a2_vals,
                colorscale=[[0, "rgba(0,0,0,1)"], [1, "rgba(0,0,0,1)"]],  # invisible
                contours=dict(start=tol, end=tol, coloring="lines"),
                showscale=False,
            )
        )
        fig.add_trace(
            go.Contour(
                z=Z_N2,
                x=a1_vals,
                y=a2_vals,
                colorscale=[[0, "rgba(0,0,0,1)"], [1, "rgba(0,0,0,1)"]],  # invisible
                contours=dict(start=tol, end=tol, coloring="lines"),
                showscale=False,
            )
        )

        fig.update_layout(
            xaxis=dict(title="a At", ticks="inside"),
            yaxis=dict(title="a Oa", ticks="inside"),
            # title=f"Steady-state {scale_names[i]}",
        )
        fig = style_plot(fig, line_thickness=1.5, font_size=font_size)
        fig.add_annotation(
            x=0.3,
            y=0.3,
            text="Oa wins",
            showarrow=False,
            xanchor="right",
            yanchor="top",
            font=dict(size=18, color="white"),
        )
        fig.add_annotation(
            x=0.8,
            y=0.8,
            text="At wins",
            showarrow=False,
            xanchor="right",
            yanchor="top",
            font=dict(size=18, color="white"),
        )

        fig.add_annotation(
            x=0.8,
            y=0.2,
            text="Coexistence",
            showarrow=False,
            xanchor="center",
            yanchor="top",
            font=dict(color="white", size=18),
        )
        fig.add_annotation(
            x=0.15,
            y=0.8,
            text="Coexistence",
            showarrow=False,
            xanchor="center",
            yanchor="top",
            font=dict(color="white", size=18),
        )
        fig.write_image(f"plots/contours/{filenames[i]}")


a_mono_culture()
a1_a2_sweep()
