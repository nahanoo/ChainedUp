from scipy.integrate import odeint
import numpy as np
import plotly.graph_objects as go
from style import *
import pandas as pd
import statsmodels.api as sm
from scipy.signal import medfilt, savgol_filter
from scipy.interpolate import UnivariateSpline


def smooth_growth_rate(signal, median_kernel=5, sg_window=11, sg_poly=3):
    signal = np.asarray(signal)

    # Step 1: median filter
    filtered = medfilt(signal, kernel_size=median_kernel)

    # Step 2: Savitzky-Golay filter
    smoothed = savgol_filter(filtered, window_length=sg_window, polyorder=sg_poly)

    return smoothed


lowess = sm.nonparametric.lowess
xs = np.linspace(0, 10, 60)
meta = pd.read_csv("data/250918_glucose_lactose_plate/metadata.csv")
data = pd.read_csv("data/250918_glucose_lactose_plate/measurements.csv", index_col=0)


def get_linegroups(species, cs, signal, meta):
    linegroups = meta[
        (meta["species"] == species)
        & (meta["carbon_source"] == cs)
        & (meta["exp_ID"] == signal)
    ]["linegroup"].unique()
    return linegroups


def model(y, t):
    N1, N2, R = y
    J1 = 0.4 * R / (0.1 + R)
    J2 = 0.2 * R / (0.01 + R)
    dN1 = J1 * N1
    dN2 = J2 * N2
    dR = -J1 * N1 / 1 - J2 * N2 / 0.5
    return [dN1, dN2, dR]


def add_noise(N):
    noise_level = 0.05
    noise = np.random.normal(0, noise_level * N.max(), size=N.shape)
    return N + noise


def simulate():
    N1, N2, R = odeint(model, [0.1, 0.1, 1], xs).T
    N1 = add_noise(N1)
    N2 = add_noise(N2)
    total = add_noise(N1 + N2)
    return N1, N2, total


def plot_simulation():
    N1, N2, total = simulate()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xs, y=N1, mode="lines", name="N1"))
    fig.add_trace(go.Scatter(x=xs, y=N2, mode="lines", name="N2"))
    fig.add_trace(go.Scatter(x=xs, y=total, mode="lines", name="total"))
    fig.update_layout(
        title="Population Dynamics with Noise",
        xaxis_title="Time",
        yaxis_title="Population",
    )
    fig = style_plot(fig)
    fig.write_image("plots/simulation.svg")
    return fig


def simulate_reconstruction():
    N1, N2, total = simulate()
    r = np.diff(np.log(N1)) / np.diff(xs)
    rN1 = r * N1[:-1]

    N_reconstructed = [N1[0]]
    for i, m in enumerate(r):
        dt = xs[i + 1] - xs[i]
        N_next = N_reconstructed[-1] * np.exp(m * dt)
        N_reconstructed.append(N_next)
    N_reconstructed = np.array(N_reconstructed)

    N2_reconstructed = total - N_reconstructed

    fig = plot_simulation()
    fig.add_trace(go.Scatter(x=xs[:-1], y=N_reconstructed, mode="lines", name="rN1"))
    fig.add_trace(go.Scatter(x=xs[:-1], y=N2_reconstructed, mode="lines", name="rN2"))
    fig.write_image("plots/simulation_with_growth_rate.svg")


def gfp_od_linearity():
    species = "At"
    cs = "Lactose"
    signal = "glucose_lactose_250918_glucose_lactose_screening_OD"
    fig = go.Figure()
    lgs_OD = get_linegroups(species, cs, signal, meta)
    for i, od in enumerate(lgs_OD):
        x = data[od + "_time"].to_numpy()
        y_od = data[od + "_measurement"].to_numpy()
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y_od,
                mode="lines",
                name=od,
                marker=dict(color="blue"),
                showlegend=False,
            )
        )
    fig.update_layout(
        title=f"OD600 of {species} on {cs}",
        xaxis_title="Time (h)",
        yaxis_title="OD600",
    )
    fig = style_plot(fig)
    fig.write_image("plots/250918_glucose_lactose_od.svg")
    lgs_OD = get_linegroups(species, cs, signal, meta)
    signal = "glucose_lactose_250918_glucose_lactose_screening_GFP"
    lgs_GFP = get_linegroups(species, cs, signal, meta)
    fig = go.Figure()
    for i, gfp in enumerate(lgs_GFP):
        x = data[gfp + "_time"].to_numpy()
        y_gfp = data[gfp + "_measurement"].to_numpy()
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y_gfp,
                mode="lines",
                name=gfp,
                marker=dict(color="red"),
                showlegend=False,
            )
        )
    fig.update_layout(
        title=f"GFP of {species} on {cs}",
        xaxis_title="Time (h)",
        yaxis_title="GFP",
    )
    fig = style_plot(fig)
    fig.write_image("plots/250918_glucose_lactose_gfp.svg")
    fig = go.Figure()
    for i, (od, gfp) in enumerate(zip(lgs_OD, lgs_GFP)):
        x = data[od + "_time"].to_numpy()
        y_od = data[od + "_measurement"].to_numpy()[1:-1]
        y_gfp = data[gfp + "_measurement"].to_numpy()[1:-1]
        keep_od = []
        keep_y_gfp = []
        keep_x = []
        for x_ods, y_ods, y_gfps in zip(x, y_od, y_gfp):
            if y_ods > 0 and y_gfps > 0:
                keep_od.append(y_ods)
                keep_y_gfp.append(y_gfps)
                keep_x.append(x_ods)
        x = np.array(keep_x)
        y_od = np.array(keep_od)
        y_gfp = np.array(keep_y_gfp)
        ratio = y_gfp / y_od
        fig.add_trace(
            go.Scatter(
                x=x,
                y=ratio,
                mode="lines",
                name=od,
                marker=dict(color="green"),
                showlegend=False,
            )
        )
    fig.update_layout(
        title=f"GFP/OD600 ratio of {species} on {cs}",
        xaxis_title="Time (h)",
        yaxis_title="GFP/OD600",
    )
    fig = style_plot(fig)
    fig.write_image("plots/250918_glucose_lactose_gfp_od_linearity.svg")


def gfp_od_inference():
    pass


species = "At"
cs = "Glucose+Lactose"
signal = "glucose_lactose_250918_glucose_lactose_screening_OD"
fig = go.Figure()
lgs_OD = get_linegroups(species, cs, signal, meta)
for i, od in enumerate(lgs_OD):
    x = data[od + "_time"].to_numpy()
    y_od = data[od + "_measurement"].to_numpy()[1:-1]
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y_od,
            mode="lines",
            name=od,
            marker=dict(color="blue"),
            showlegend=False,
        )
    )
signal = "glucose_lactose_250918_glucose_lactose_screening_GFP"
lgs_GFP = get_linegroups(species, cs, signal, meta)
fig2 = go.Figure()
fig3 = go.Figure()
i = 0
y_gfps = []
for i, gfp in enumerate(lgs_GFP):
    x = data[gfp + "_time"].to_numpy()[0:]
    y_gfp = data[gfp + "_measurement"].to_numpy()[0:]
    fig3.add_trace(
        go.Scatter(
            x=x,
            y=y_gfp,
            mode="lines",
            name=gfp,
            marker=dict(color=colors["at"]),
            showlegend=False,
        )
    )
    y_gfps.append(y_gfp)


y_gfp = np.mean(y_gfps, axis=0)
y_gfp_keep = []
x_keep = []
for i, (xi, yi) in enumerate(zip(x, y_gfp)):
    if yi > 0:
        y_gfp_keep.append(yi)
        x_keep.append(xi)
y_gfp = np.array(y_gfp_keep)
x = np.array(x_keep)
r_gfp = np.diff(np.log(y_gfp)) / np.diff(x) / 1.4
fig2.add_trace(
    go.Scatter(
        x=x[:-1],
        y=r_gfp,
        mode="lines",
        name=gfp[-2:],
        # marker=dict(color="red"),
        showlegend=True,
    )
)
N_reconstructed = [y_od[0]]
for i, m in enumerate(r_gfp[:-1]):
    dt = x[i + 1] - x[i]
    N_next = N_reconstructed[-1] * np.exp(m * dt)
    N_reconstructed.append(N_next)
N_reconstructed = np.array(N_reconstructed)
fig.add_trace(
    go.Scatter(
        x=x,
        y=N_reconstructed,
        mode="lines",
        name=gfp,
        marker=dict(color=colors["at"]),
        showlegend=False,
    )
)
fig2.update_layout(
    title=f"GFP inferred growth rate of {species} on {cs}",
    xaxis_title="Time (h)",
    yaxis_title="Inferred growth rate (1/h)",
)
fig2 = style_plot(fig2)
fig2.write_image("plots/250918_glucose_lactose_gfp_growth_rate.svg")

fig3.update_layout(xaxis_title="Time (h)", yaxis_title="GFP")
fig3 = style_plot(fig3)
fig3.write_image("plots/250918_glucose_lactose_gfp.svg")
fig.update_layout(
    title=f"GFP inferred growth of {species} on {cs}",
    xaxis_title="Time (h)",
    yaxis_title="Inferred OD600",
)
fig = style_plot(fig)
fig.write_image("plots/250918_glucose_lactose_gfp_inference.svg")

fig = gfp_od_inference()

species = "At"
cs = "Lactose"
signal = "glucose_lactose_250918_glucose_lactose_screening_OD"
fig = go.Figure()
lgs_OD = get_linegroups(species, cs, signal, meta)
