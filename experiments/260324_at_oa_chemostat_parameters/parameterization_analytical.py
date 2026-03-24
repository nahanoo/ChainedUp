import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.integrate import solve_ivp

import conditions as cond
from style import *


DILUTION_RATE = 0.15
N_STEADY_STATE_POINTS = 3

OUTPUT_DIR = "plots/parameterization"
PARAMETER_CSV = "parameterization_results_analytical_peak_area_units.csv"


SPECS = {
    "succinate": {
        "wells": cond.OA_SUCCINATE,
        "times": np.array(cond.OA_TIMEPOINTS, dtype=float),
        "substrate_col": "succinate_area_normalized",
        "reactor": "Succinate",
        "mumax": 0.877,
        "substrate_label": "Succinate",
        "substrate_color_key": "Succinate",
        "fname": "oa_succinate_analytical_fit.svg",
    },
    "glucose": {
        "wells": cond.OA_GLUCOSE,
        "times": np.array(cond.OA_TIMEPOINTS, dtype=float),
        "substrate_col": "glucose_area_normalized",
        "reactor": "Glucose",
        "mumax": 0.573,
        "substrate_label": "Glucose",
        "substrate_color_key": "Glucose",
        "fname": "oa_glucose_analytical_fit.svg",
    },
}


def monod_chemostat(t, y, mumax, km, yld, dilution_rate, s_in):
    n, s = y

    mu = mumax * s / (km + s)

    dn_dt = (mu - dilution_rate) * n
    ds_dt = dilution_rate * (s_in - s) - (1 / yld) * mu * n

    return [dn_dt, ds_dt]


def simulate(times, n0, s0, mumax, km, yld, dilution_rate, s_in):
    sol = solve_ivp(
        fun=lambda t, y: monod_chemostat(
            t=t,
            y=y,
            mumax=mumax,
            km=km,
            yld=yld,
            dilution_rate=dilution_rate,
            s_in=s_in,
        ),
        t_span=(times[0], times[-1]),
        y0=[n0, s0],
        t_eval=times,
        method="LSODA",
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    return sol.y[0], sol.y[1]


def prepare_condition_data(resource_df, cfu_df, wells, times, substrate_col, reactor):
    substrate_df = pd.DataFrame(
        {
            "time": times,
            "s_obs": [resource_df.loc[well, substrate_col] for well in wells],
        }
    )

    biomass_df = (
        cfu_df[(cfu_df["species"].str.lower() == "oa") & (cfu_df["reactor"] == reactor)]
        .groupby("sample_time", as_index=False)["average"]
        .mean()
        .rename(columns={"sample_time": "time", "average": "n_obs"})
    )

    merged = (
        biomass_df.merge(substrate_df, on="time", how="inner")
        .sort_values("time")
        .reset_index(drop=True)
    )

    if merged.empty:
        raise ValueError(f"No overlapping timepoints found for reactor {reactor!r}.")

    return (
        merged["time"].to_numpy(dtype=float),
        merged["n_obs"].to_numpy(dtype=float),
        merged["s_obs"].to_numpy(dtype=float),
    )


def compute_analytical_parameters(
    times,
    n_obs,
    s_obs,
    mumax,
    dilution_rate,
    n_steady_state_points=N_STEADY_STATE_POINTS,
):
    if len(times) < n_steady_state_points:
        raise ValueError(
            f"Need at least {n_steady_state_points} timepoints, got {len(times)}."
        )

    if mumax <= dilution_rate:
        raise ValueError(
            f"mumax ({mumax}) must be larger than dilution rate ({dilution_rate})."
        )

    steady_idx = np.arange(len(times) - n_steady_state_points, len(times))

    n_star = float(np.mean(n_obs[steady_idx]))
    s_star = float(np.mean(s_obs[steady_idx]))

    n_star_std = (
        float(np.std(n_obs[steady_idx], ddof=1)) if len(steady_idx) > 1 else 0.0
    )
    s_star_std = (
        float(np.std(s_obs[steady_idx], ddof=1)) if len(steady_idx) > 1 else 0.0
    )

    s_in = float(s_obs[0])

    km = s_star * (mumax - dilution_rate) / dilution_rate

    delta_s = s_in - s_star
    if delta_s <= 0:
        raise ValueError(
            f"s_in - s_star must be positive to calculate yield, got {delta_s:.6g}."
        )

    yld = n_star / delta_s

    return {
        "steady_idx": steady_idx,
        "n_star": n_star,
        "s_star": s_star,
        "n_star_std": n_star_std,
        "s_star_std": s_star_std,
        "s_in": s_in,
        "km": km,
        "yield": yld,
    }


def make_analytical_plot(name, fit, spec):
    times = fit["times"]

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=("Biomass", spec["substrate_label"]),
    )

    biomass_color = colors["Oa"] if "Oa" in colors else qualitative.Safe[0]
    substrate_color = colors[spec["substrate_color_key"]]

    fig.add_trace(
        go.Scatter(
            x=times,
            y=fit["n_obs"],
            mode="markers",
            name="CFU data",
            marker=dict(
                color=biomass_color,
                size=8,
                line=dict(color="black", width=1.2),
            ),
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=fit["t_sim"],
            y=fit["n_sim"],
            mode="lines",
            name="CFU simulation",
            line=dict(color=biomass_color, width=3),
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=times[fit["steady_idx"]],
            y=fit["n_obs"][fit["steady_idx"]],
            mode="markers",
            name="CFU steady-state points",
            marker=dict(
                color=biomass_color,
                size=11,
                symbol="circle-open",
                line=dict(color="black", width=2),
            ),
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=times,
            y=fit["s_obs"],
            mode="markers",
            name=f"{spec['substrate_label']} data",
            marker=dict(
                color=substrate_color,
                size=8,
                line=dict(color="black", width=1.2),
            ),
        ),
        row=2,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=fit["t_sim"],
            y=fit["s_sim"],
            mode="lines",
            name=f"{spec['substrate_label']} simulation",
            line=dict(color=substrate_color, width=3),
        ),
        row=2,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=times[fit["steady_idx"]],
            y=fit["s_obs"][fit["steady_idx"]],
            mode="markers",
            name=f"{spec['substrate_label']} steady-state points",
            marker=dict(
                color=substrate_color,
                size=11,
                symbol="circle-open",
                line=dict(color="black", width=2),
            ),
        ),
        row=2,
        col=1,
    )

    fig.add_hline(
        y=fit["s_in"],
        line=dict(color=substrate_color, width=2, dash="dot"),
        row=2,
        col=1,
    )

    fig.update_yaxes(
        title_text="CFU/mL",
        type="log",
        exponentformat="power",
        row=1,
        col=1,
    )
    fig.update_yaxes(
        title_text="Normalized peak area",
        type="log",
        exponentformat="power",
        row=2,
        col=1,
    )
    fig.update_xaxes(title_text="Time (h)", row=2, col=1)

    fig.update_layout(
        width=380,
        height=430,
        title=(
            f"Oa {name} analytical parameterization"
            f"<br><sup>"
            f"Km={fit['km']:.3g} peak area, "
            f"Yield={fit['yield']:.3g} CFU/peak area, "
            f"mu_max={spec['mumax']}, "
            f"D={DILUTION_RATE}"
            f"</sup>"
        ),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255,255,255,0)",
        ),
    )

    fig = style_plot(fig, marker_size=8, line_thickness=3, font_size=11)
    return fig


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    resource_df = pd.read_csv("normalized_data.csv", index_col="id")
    cfu_df = pd.read_csv("cfus.csv")

    results = []

    for name, spec in SPECS.items():
        times, n_obs, s_obs = prepare_condition_data(
            resource_df=resource_df,
            cfu_df=cfu_df,
            wells=spec["wells"],
            times=spec["times"],
            substrate_col=spec["substrate_col"],
            reactor=spec["reactor"],
        )

        analytical = compute_analytical_parameters(
            times=times,
            n_obs=n_obs,
            s_obs=s_obs,
            mumax=spec["mumax"],
            dilution_rate=DILUTION_RATE,
            n_steady_state_points=N_STEADY_STATE_POINTS,
        )

        n0 = float(n_obs[0])
        s0 = float(s_obs[0])

        t_sim = np.linspace(times[0], times[-1], 400)
        n_sim, s_sim = simulate(
            times=t_sim,
            n0=n0,
            s0=s0,
            mumax=spec["mumax"],
            km=analytical["km"],
            yld=analytical["yield"],
            dilution_rate=DILUTION_RATE,
            s_in=analytical["s_in"],
        )

        fit = {
            "times": times,
            "n_obs": n_obs,
            "s_obs": s_obs,
            "steady_idx": analytical["steady_idx"],
            "n_star": analytical["n_star"],
            "s_star": analytical["s_star"],
            "n_star_std": analytical["n_star_std"],
            "s_star_std": analytical["s_star_std"],
            "s_in": analytical["s_in"],
            "km": analytical["km"],
            "yield": analytical["yield"],
            "n0": n0,
            "s0": s0,
            "t_sim": t_sim,
            "n_sim": n_sim,
            "s_sim": s_sim,
        }

        results.append(
            {
                "condition": name,
                "mumax_h-1": spec["mumax"],
                "dilution_rate_h-1": DILUTION_RATE,
                "n_steady_state_points": N_STEADY_STATE_POINTS,
                "Km_peak_area": fit["km"],
                "yield_cfu_per_peak_area": fit["yield"],
                "n0_cfu_ml": fit["n0"],
                "s0_peak_area": fit["s0"],
                "sin_peak_area": fit["s_in"],
                "n_star_cfu_ml": fit["n_star"],
                "n_star_std_cfu_ml": fit["n_star_std"],
                "s_star_peak_area": fit["s_star"],
                "s_star_std_peak_area": fit["s_star_std"],
            }
        )

        fig = make_analytical_plot(name=name, fit=fit, spec=spec)
        fig.write_image(os.path.join(OUTPUT_DIR, spec["fname"]))

    results_df = pd.DataFrame(results)
    results_df.to_csv(PARAMETER_CSV, index=False)

    print(results_df.to_string(index=False))


if __name__ == "__main__":
    main()
