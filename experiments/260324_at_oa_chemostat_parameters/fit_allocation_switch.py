import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import conditions as cond
from style import *


DILUTION_RATE = 0.3
TIME_MAX = 100
SWITCH_TIME = 52.8

MONO_PARAMETER_CSV = "parameterization_results.csv"
OUTPUT_DIR = "plots/allocation"
RESULT_CSV = "allocation_switch_fit_result.csv"


MIXED_SPEC = {
    "wells": cond.OA_SUCCINATE_GLUCOSE,
    "times": np.array(cond.OA_TIMEPOINTS, dtype=float),
    "reactor": "Succinate+Glucose",
    "succinate_col": "succinate_area_normalized",
    "glucose_col": "glucose_area_normalized",
    "cfu_fname": "oa_succinate_glucose_0_60_cfu_switch_fit.svg",
    "resource_fname": "oa_succinate_glucose_0_60_resource_switch_fit.svg",
}


def load_mono_parameters(parameter_csv):
    params_df = pd.read_csv(parameter_csv)
    params_df["condition"] = params_df["condition"].str.lower()

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


def simulate_segment(
    t_start, y0, times_eval, a, mono_params, dilution_rate, s_in, g_in
):
    if len(times_eval) == 0:
        return (
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
        )

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
        t_span=(t_start, times_eval[-1]),
        y0=y0,
        t_eval=times_eval,
        method="LSODA",
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    return sol.y[0], sol.y[1], sol.y[2]


def simulate_state_at_time(
    t_start, y0, t_end, a, mono_params, dilution_rate, s_in, g_in
):
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
        t_span=(t_start, t_end),
        y0=y0,
        t_eval=[t_end],
        method="LSODA",
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    return sol.y[:, -1]


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
        raise ValueError(
            "No overlapping mixed-condition timepoints found up to TIME_MAX."
        )

    return {
        "times": merged["time"].to_numpy(dtype=float),
        "n_obs": merged["n_obs"].to_numpy(dtype=float),
        "succinate_obs": merged["succinate_obs"].to_numpy(dtype=float),
        "glucose_obs": merged["glucose_obs"].to_numpy(dtype=float),
    }


def residuals_segment(
    a_array,
    t_start,
    y0,
    times_eval,
    n_obs,
    s_obs,
    g_obs,
    mono_params,
    dilution_rate,
    s_in,
    g_in,
):
    a = float(a_array[0])

    n_pred, s_pred, g_pred = simulate_segment(
        t_start=t_start,
        y0=y0,
        times_eval=times_eval,
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


def fit_segment(
    t_start,
    y0,
    times_eval,
    n_obs,
    s_obs,
    g_obs,
    mono_params,
    dilution_rate,
    s_in,
    g_in,
):
    result = least_squares(
        fun=residuals_segment,
        x0=np.array([0.5]),
        bounds=(np.array([1e-6]), np.array([1 - 1e-6])),
        args=(
            t_start,
            y0,
            times_eval,
            n_obs,
            s_obs,
            g_obs,
            mono_params,
            dilution_rate,
            s_in,
            g_in,
        ),
    )

    a = float(result.x[0])

    n_pred, s_pred, g_pred = simulate_segment(
        t_start=t_start,
        y0=y0,
        times_eval=times_eval,
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
        "n_pred": n_pred,
        "s_pred": s_pred,
        "g_pred": g_pred,
    }


def fit_two_phase_allocation(mixed_data, mono_params, dilution_rate, switch_time):
    times = mixed_data["times"]
    n_obs = mixed_data["n_obs"]
    s_obs = mixed_data["succinate_obs"]
    g_obs = mixed_data["glucose_obs"]

    pre_mask = times <= switch_time
    post_mask = times > switch_time

    if pre_mask.sum() < 2:
        raise ValueError("Not enough pre-switch timepoints to fit a_pre.")
    if post_mask.sum() < 1:
        raise ValueError("No post-switch timepoints available to fit a_post.")

    times_pre = times[pre_mask]
    n_pre = n_obs[pre_mask]
    s_pre = s_obs[pre_mask]
    g_pre = g_obs[pre_mask]

    times_post = times[post_mask]
    n_post = n_obs[post_mask]
    s_post = s_obs[post_mask]
    g_post = g_obs[post_mask]

    n0 = float(n_obs[0])
    s0 = float(s_obs[0])
    g0 = float(g_obs[0])

    y0 = [n0, s0, g0]
    s_in = s0
    g_in = g0

    pre_fit = fit_segment(
        t_start=times_pre[0],
        y0=y0,
        times_eval=times_pre,
        n_obs=n_pre,
        s_obs=s_pre,
        g_obs=g_pre,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    state_switch = simulate_state_at_time(
        t_start=times_pre[0],
        y0=y0,
        t_end=switch_time,
        a=pre_fit["a"],
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    post_fit = fit_segment(
        t_start=switch_time,
        y0=state_switch,
        times_eval=times_post,
        n_obs=n_post,
        s_obs=s_post,
        g_obs=g_post,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    t_pre_dense = np.linspace(times[0], switch_time, 250)
    t_post_dense = np.linspace(switch_time, TIME_MAX, 250)

    n_pre_sim, s_pre_sim, g_pre_sim = simulate_segment(
        t_start=times[0],
        y0=y0,
        times_eval=t_pre_dense,
        a=pre_fit["a"],
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    state_switch_dense = np.array([n_pre_sim[-1], s_pre_sim[-1], g_pre_sim[-1]])

    n_post_sim, s_post_sim, g_post_sim = simulate_segment(
        t_start=switch_time,
        y0=state_switch_dense,
        times_eval=t_post_dense,
        a=post_fit["a"],
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    t_sim = np.concatenate([t_pre_dense, t_post_dense[1:]])
    n_sim = np.concatenate([n_pre_sim, n_post_sim[1:]])
    s_sim = np.concatenate([s_pre_sim, s_post_sim[1:]])
    g_sim = np.concatenate([g_pre_sim, g_post_sim[1:]])

    return {
        "times": times,
        "n_obs": n_obs,
        "succinate_obs": s_obs,
        "glucose_obs": g_obs,
        "a_pre": pre_fit["a"],
        "a_post": post_fit["a"],
        "pre_success": pre_fit["success"],
        "post_success": post_fit["success"],
        "pre_message": pre_fit["message"],
        "post_message": post_fit["message"],
        "pre_cost": pre_fit["cost"],
        "post_cost": post_fit["cost"],
        "switch_time": switch_time,
        "n0": n0,
        "s0": s0,
        "g0": g0,
        "s_in": s_in,
        "g_in": g_in,
        "t_sim": t_sim,
        "n_sim": n_sim,
        "succinate_sim": s_sim,
        "glucose_sim": g_sim,
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

    fig.add_vline(
        x=fit["switch_time"],
        line=dict(color="gray", dash="dash", width=1.5),
    )

    fig.update_layout(
        width=200,
        height=180,
        title="Oa Succinate Glucose",
        showlegend=False,
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

    fig.add_vline(
        x=fit["switch_time"],
        line=dict(color="gray", dash="dash", width=1.5),
    )

    fig.update_layout(
        height=200,
        width=260,
        title="Oa Succinate Glucose",
        showlegend=False,
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

    fit = fit_two_phase_allocation(
        mixed_data=mixed_data,
        mono_params=mono_params,
        dilution_rate=DILUTION_RATE,
        switch_time=SWITCH_TIME,
    )

    result_row = {
        "switch_time": fit["switch_time"],
        "time_max": TIME_MAX,
        "dilution_rate": DILUTION_RATE,
        "a_pre_succinate": fit["a_pre"],
        "a_pre_glucose": 1 - fit["a_pre"],
        "a_post_succinate": fit["a_post"],
        "a_post_glucose": 1 - fit["a_post"],
        "pre_success": fit["pre_success"],
        "post_success": fit["post_success"],
        "pre_cost": fit["pre_cost"],
        "post_cost": fit["post_cost"],
        "mumax_succinate": mono_params["succinate"]["mumax"],
        "km_succinate": mono_params["succinate"]["km"],
        "yield_succinate": mono_params["succinate"]["yield"],
        "mumax_glucose": mono_params["glucose"]["mumax"],
        "km_glucose": mono_params["glucose"]["km"],
        "yield_glucose": mono_params["glucose"]["yield"],
        "pre_message": fit["pre_message"],
        "post_message": fit["post_message"],
    }

    results_df = pd.DataFrame([result_row])
    results_df.to_csv(RESULT_CSV, index=False)

    fig_cfu = make_cfu_plot(fit)
    fig_cfu.write_image(os.path.join(OUTPUT_DIR, MIXED_SPEC["cfu_fname"]))

    fig_resource = make_resource_plot(fit)
    fig_resource.write_image(os.path.join(OUTPUT_DIR, MIXED_SPEC["resource_fname"]))

    print(results_df.to_string(index=False))
    print(f"\na_pre  (succinate allocation) = {fit['a_pre']:.6f}")
    print(f"a_pre  (glucose allocation)   = {1 - fit['a_pre']:.6f}")
    print(f"a_post (succinate allocation) = {fit['a_post']:.6f}")
    print(f"a_post (glucose allocation)   = {1 - fit['a_post']:.6f}")


if __name__ == "__main__":
    main()
