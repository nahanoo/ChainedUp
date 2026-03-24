import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import conditions as cond
from style import *


DILUTION_RATE = 0.3
TIME_MAX = 60

MONO_PARAMETER_CSV = "parameterization_results.csv"
OUTPUT_DIR = "plots/allocation"
RESULT_CSV = "allocation_fit_result.csv"


MIXED_SPEC = {
    "wells": cond.OA_SUCCINATE_GLUCOSE,
    "times": np.array(cond.OA_TIMEPOINTS, dtype=float),
    "reactor": "Succinate+Glucose",
    "succinate_col": "succinate_area_normalized",
    "glucose_col": "glucose_area_normalized",
    "cfu_fname": "oa_succinate_glucose_0_60_cfu_allocation_fit.svg",
    "resource_fname": "oa_succinate_glucose_0_60_resource_allocation_fit.svg",
}


def load_mono_parameters(parameter_csv):
    params_df = pd.read_csv(parameter_csv)
    params_df["condition"] = params_df["condition"].str.lower()

    required_conditions = {"succinate", "glucose"}
    found_conditions = set(params_df["condition"])
    missing = required_conditions - found_conditions
    if missing:
        raise ValueError(f"Missing mono-substrate fits for: {sorted(missing)}")

    mono_params = {}
    for condition in ["succinate", "glucose"]:
        row = params_df.loc[params_df["condition"] == condition].iloc[0]
        mono_params[condition] = {
            "mumax": float(row["mumax"]),
            "km": float(row["Km"]),
            "yield": float(row["yield"]),
        }

    return mono_params


def allocation_chemostat(t, y, a, mono_params, dilution_rate, s_in, g_in):
    n, s, g = y

    mu_s = (
        a * mono_params["succinate"]["mumax"] * s / (mono_params["succinate"]["km"] + s)
    )
    mu_g = (
        (1 - a)
        * mono_params["glucose"]["mumax"]
        * g
        / (mono_params["glucose"]["km"] + g)
    )

    dn_dt = (mu_s + mu_g - dilution_rate) * n
    ds_dt = (
        dilution_rate * (s_in - s) - (1 / mono_params["succinate"]["yield"]) * mu_s * n
    )
    dg_dt = (
        dilution_rate * (g_in - g) - (1 / mono_params["glucose"]["yield"]) * mu_g * n
    )

    return [dn_dt, ds_dt, dg_dt]


def simulate_mixed(times, n0, s0, g0, a, mono_params, dilution_rate, s_in, g_in):
    sol = solve_ivp(
        fun=lambda t, y: allocation_chemostat(
            t=t,
            y=y,
            a=a,
            mono_params=mono_params,
            dilution_rate=dilution_rate,
            s_in=s_in,
            g_in=g_in,
        ),
        t_span=(times[0], times[-1]),
        y0=[n0, s0, g0],
        t_eval=times,
        method="LSODA",
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    return sol.y[0], sol.y[1], sol.y[2]


def prepare_mixed_data(resource_df, cfu_df):
    resource_sub = pd.DataFrame(
        {
            "time": MIXED_SPEC["times"],
            "succinate_obs": [
                resource_df.loc[well, MIXED_SPEC["succinate_col"]]
                for well in MIXED_SPEC["wells"]
            ],
            "glucose_obs": [
                resource_df.loc[well, MIXED_SPEC["glucose_col"]]
                for well in MIXED_SPEC["wells"]
            ],
        }
    )

    biomass_df = (
        cfu_df[
            (cfu_df["species"].str.lower() == "oa")
            & (cfu_df["reactor"] == MIXED_SPEC["reactor"])
        ]
        .groupby("sample_time", as_index=False)["average"]
        .mean()
        .rename(columns={"sample_time": "time", "average": "n_obs"})
    )

    merged = (
        biomass_df.merge(resource_sub, on="time", how="inner")
        .sort_values("time")
        .reset_index(drop=True)
    )

    merged = merged[merged["time"] <= TIME_MAX].reset_index(drop=True)

    if merged.empty:
        raise ValueError("No overlapping mixed-condition timepoints found up to 60 h.")

    return {
        "times": merged["time"].to_numpy(dtype=float),
        "n_obs": merged["n_obs"].to_numpy(dtype=float),
        "succinate_obs": merged["succinate_obs"].to_numpy(dtype=float),
        "glucose_obs": merged["glucose_obs"].to_numpy(dtype=float),
    }


def residuals(
    a_array, times, n_obs, s_obs, g_obs, mono_params, dilution_rate, s_in, g_in
):
    a = float(a_array[0])

    n0 = float(n_obs[0])
    s0 = float(s_obs[0])
    g0 = float(g_obs[0])

    n_pred, s_pred, g_pred = simulate_mixed(
        times=times,
        n0=n0,
        s0=s0,
        g0=g0,
        a=a,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    eps = 1e-12

    n_res = np.log10(np.clip(n_pred, eps, None)) - np.log10(np.clip(n_obs, eps, None))
    s_res = np.log10(np.clip(s_pred, eps, None)) - np.log10(np.clip(s_obs, eps, None))
    g_res = np.log10(np.clip(g_pred, eps, None)) - np.log10(np.clip(g_obs, eps, None))

    return np.concatenate([n_res, s_res, g_res])


def fit_allocation_parameter(times, n_obs, s_obs, g_obs, mono_params, dilution_rate):
    n0 = float(n_obs[0])
    s0 = float(s_obs[0])
    g0 = float(g_obs[0])

    s_in = s0
    g_in = g0

    result = least_squares(
        fun=residuals,
        x0=np.array([0.5]),
        bounds=(np.array([1e-6]), np.array([1 - 1e-6])),
        args=(times, n_obs, s_obs, g_obs, mono_params, dilution_rate, s_in, g_in),
    )

    a = float(result.x[0])

    n_pred, s_pred, g_pred = simulate_mixed(
        times=times,
        n0=n0,
        s0=s0,
        g0=g0,
        a=a,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    t_sim = np.linspace(times[0], times[-1], 400)
    n_sim, s_sim, g_sim = simulate_mixed(
        times=t_sim,
        n0=n0,
        s0=s0,
        g0=g0,
        a=a,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    return {
        "a": a,
        "success": bool(result.success),
        "message": result.message,
        "cost": float(result.cost),
        "times": times,
        "n_obs": n_obs,
        "succinate_obs": s_obs,
        "glucose_obs": g_obs,
        "n_pred": n_pred,
        "succinate_pred": s_pred,
        "glucose_pred": g_pred,
        "t_sim": t_sim,
        "n_sim": n_sim,
        "succinate_sim": s_sim,
        "glucose_sim": g_sim,
        "n0": n0,
        "s0": s0,
        "g0": g0,
        "s_in": s_in,
        "g_in": g_in,
    }


def make_cfu_plot(fit):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["n_obs"],
            mode="markers",
            name="Data",
            marker=dict(
                color=colors["Succinate+Glucose"],
                size=8,
                opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=fit["t_sim"],
            y=fit["n_sim"],
            mode="lines",
            name="Fit",
            line=dict(color=colors["Succinate+Glucose"], width=3),
            opacity=0.8,
        )
    )

    fig.update_layout(
        width=200,
        height=180,
        title="Oa Succinate Glucose",
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95,
            bgcolor="rgba(255,255,255,0)",
        ),
        xaxis=dict(title="Time (h)", range=[0, TIME_MAX], dtick=20),
        yaxis=dict(
            title="CFU/mL",
            type="log",
            exponentformat="power",
            range=[6, 10],
        ),
    )

    fig = style_plot(fig, marker_size=8, line_thickness=3, font_size=11)
    return fig


def make_resource_plot(fit):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["succinate_obs"],
            mode="markers",
            name="Succinate data",
            marker=dict(
                color=colors["Succinate"],
                size=8,
                opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=fit["t_sim"],
            y=fit["succinate_sim"],
            mode="lines",
            name="Succinate fit",
            line=dict(color=colors["Succinate"], width=3),
            opacity=0.8,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["glucose_obs"],
            mode="markers",
            name="Glucose data",
            marker=dict(
                color=colors["Glucose"],
                size=8,
                opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=fit["t_sim"],
            y=fit["glucose_sim"],
            mode="lines",
            name="Glucose fit",
            line=dict(color=colors["Glucose"], width=3),
            opacity=0.8,
        )
    )

    fig.update_layout(
        width=200,
        height=180,
        title="Oa Succinate Glucose",
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95,
            bgcolor="rgba(255,255,255,0)",
        ),
        xaxis=dict(title="Time (h)", range=[0, TIME_MAX], dtick=20),
        yaxis=dict(
            title="Area",
            type="log",
            exponentformat="power",
        ),
    )

    fig = style_plot(fig, marker_size=8, line_thickness=3, font_size=11)
    return fig


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    mono_params = load_mono_parameters(MONO_PARAMETER_CSV)
    resource_df = pd.read_csv("normalized_data.csv", index_col="id")
    cfu_df = pd.read_csv("cfus.csv")

    mixed_data = prepare_mixed_data(resource_df=resource_df, cfu_df=cfu_df)

    fit = fit_allocation_parameter(
        times=mixed_data["times"],
        n_obs=mixed_data["n_obs"],
        s_obs=mixed_data["succinate_obs"],
        g_obs=mixed_data["glucose_obs"],
        mono_params=mono_params,
        dilution_rate=DILUTION_RATE,
    )

    result_row = {
        "a_succinate": fit["a"],
        "a_glucose": 1 - fit["a"],
        "dilution_rate": DILUTION_RATE,
        "time_max": TIME_MAX,
        "success": fit["success"],
        "cost": fit["cost"],
        "message": fit["message"],
        "mumax_succinate": mono_params["succinate"]["mumax"],
        "km_succinate": mono_params["succinate"]["km"],
        "yield_succinate": mono_params["succinate"]["yield"],
        "mumax_glucose": mono_params["glucose"]["mumax"],
        "km_glucose": mono_params["glucose"]["km"],
        "yield_glucose": mono_params["glucose"]["yield"],
    }

    results_df = pd.DataFrame([result_row])
    results_df.to_csv(RESULT_CSV, index=False)

    fig_cfu = make_cfu_plot(fit)
    fig_cfu.write_image(os.path.join(OUTPUT_DIR, MIXED_SPEC["cfu_fname"]))

    fig_resource = make_resource_plot(fit)
    fig_resource.write_image(os.path.join(OUTPUT_DIR, MIXED_SPEC["resource_fname"]))

    print(results_df.to_string(index=False))
    print(f"\nFitted allocation parameter a = {fit['a']:.6f}")
    print(f"This means allocation to succinate = {fit['a']:.6f}")
    print(f"And allocation to glucose = {1 - fit['a']:.6f}")


if __name__ == "__main__":
    main()
