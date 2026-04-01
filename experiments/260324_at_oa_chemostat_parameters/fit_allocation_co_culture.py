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

MANUAL_A = True
A_AT = 0.5
A_OA = 0.3

# Set to False to fit one constant allocation for each species.
# Set to True to fit one allocation before SWITCH_TIME and one after.
USE_SWITCH = False
SWITCH_TIME = 40

MONO_PARAMETER_CSV = "parameterization_results.csv"
OUTPUT_DIR = "plots/allocation"


MIXED_SPEC = {
    "wells": cond.ATOA_SUCCINATE_GLUCOSE,
    "times": np.array(cond.AT_OA_TIMEPOINTS, dtype=float),
    "reactor": "Succinate+Glucose",
    "succinate_col": "succinate_area_normalized",
    "glucose_col": "glucose_area_normalized",
}


def get_output_paths(use_switch):
    mode = "switch" if use_switch else "single"
    return {
        "result_csv": f"allocation_coculture_{mode}_fit_result.csv",
        "cfu_fname": f"atoa_succinate_glucose_cfu_{mode}_fit.svg",
        "resource_fname": f"atoa_succinate_glucose_resource_{mode}_fit.svg",
    }


def load_mono_parameters(parameter_csv):
    params_df = pd.read_csv(parameter_csv)
    params_df["condition"] = params_df["condition"].str.lower()
    params_df["species"] = params_df["species"].str.lower()

    mono_params = {
        "oa": {"succinate": {}, "glucose": {}},
        "at": {"succinate": {}, "glucose": {}},
    }

    for species in ["oa", "at"]:
        for condition in ["succinate", "glucose"]:
            subset = params_df.loc[
                (params_df["species"] == species)
                & (params_df["condition"] == condition)
            ]

            if subset.empty:
                raise ValueError(
                    f"Missing row for species={species!r}, condition={condition!r} "
                    f"in {parameter_csv}."
                )
            if len(subset) > 1:
                raise ValueError(
                    f"Multiple rows found for species={species!r}, "
                    f"condition={condition!r} in {parameter_csv}."
                )

            row = subset.iloc[0]
            mono_params[species][condition] = {
                "mumax": float(row["mumax"]),
                "km": float(row["Km"]),
                "yield": float(row["yield"]),
            }

    return mono_params


def coculture_chemostat(
    t,
    y,
    a_at,
    a_oa,
    mono_params,
    dilution_rate,
    s_in,
    g_in,
):
    n_at, n_oa, s, g = y

    mu_at_s = (
        a_at
        * mono_params["at"]["succinate"]["mumax"]
        * s
        / (mono_params["at"]["succinate"]["km"] + s)
    )
    mu_at_g = (
        (1 - a_at)
        * mono_params["at"]["glucose"]["mumax"]
        * g
        / (mono_params["at"]["glucose"]["km"] + g)
    )

    mu_oa_s = (
        a_oa
        * mono_params["oa"]["succinate"]["mumax"]
        * s
        / (mono_params["oa"]["succinate"]["km"] + s)
    )
    mu_oa_g = (
        (1 - a_oa)
        * mono_params["oa"]["glucose"]["mumax"]
        * g
        / (mono_params["oa"]["glucose"]["km"] + g)
    )

    dn_at_dt = (mu_at_s + mu_at_g - dilution_rate) * n_at
    dn_oa_dt = (mu_oa_s + mu_oa_g - dilution_rate) * n_oa

    ds_dt = (
        dilution_rate * (s_in - s)
        - (1 / mono_params["at"]["succinate"]["yield"]) * mu_at_s * n_at
        - (1 / mono_params["oa"]["succinate"]["yield"]) * mu_oa_s * n_oa
    )

    dg_dt = (
        dilution_rate * (g_in - g)
        - (1 / mono_params["at"]["glucose"]["yield"]) * mu_at_g * n_at
        - (1 / mono_params["oa"]["glucose"]["yield"]) * mu_oa_g * n_oa
    )

    return [dn_at_dt, dn_oa_dt, ds_dt, dg_dt]


def simulate_segment(
    t_start,
    y0,
    times_eval,
    a_at,
    a_oa,
    mono_params,
    dilution_rate,
    s_in,
    g_in,
):
    if len(times_eval) == 0:
        return (
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
        )

    sol = solve_ivp(
        fun=lambda t, y: coculture_chemostat(
            t=t,
            y=y,
            a_at=a_at,
            a_oa=a_oa,
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

    return sol.y[0], sol.y[1], sol.y[2], sol.y[3]


def simulate_state_at_time(
    t_start,
    y0,
    t_end,
    a_at,
    a_oa,
    mono_params,
    dilution_rate,
    s_in,
    g_in,
):
    sol = solve_ivp(
        fun=lambda t, y: coculture_chemostat(
            t=t,
            y=y,
            a_at=a_at,
            a_oa=a_oa,
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


def prepare_coculture_data(resource_df, cfu_df):
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

    at_df = (
        cfu_df[
            (cfu_df["species"].str.lower() == "at_coculture")
            & (cfu_df["reactor"] == MIXED_SPEC["reactor"])
        ]
        .groupby("sample_time", as_index=False)["average"]
        .mean()
        .rename(columns={"sample_time": "time", "average": "at_obs"})
    )

    oa_df = (
        cfu_df[
            (cfu_df["species"].str.lower() == "oa_coculture")
            & (cfu_df["reactor"] == MIXED_SPEC["reactor"])
        ]
        .groupby("sample_time", as_index=False)["average"]
        .mean()
        .rename(columns={"sample_time": "time", "average": "oa_obs"})
    )

    merged = (
        resource_sub.merge(at_df, on="time", how="inner")
        .merge(oa_df, on="time", how="inner")
        .sort_values("time")
        .reset_index(drop=True)
    )

    merged = merged[merged["time"] <= TIME_MAX].reset_index(drop=True)

    if merged.empty:
        raise ValueError("No overlapping coculture timepoints found up to TIME_MAX.")

    return {
        "times": merged["time"].to_numpy(dtype=float),
        "at_obs": merged["at_obs"].to_numpy(dtype=float),
        "oa_obs": merged["oa_obs"].to_numpy(dtype=float),
        "succinate_obs": merged["succinate_obs"].to_numpy(dtype=float),
        "glucose_obs": merged["glucose_obs"].to_numpy(dtype=float),
    }


def residuals_segment(
    params_array,
    t_start,
    y0,
    times_eval,
    at_obs,
    oa_obs,
    s_obs,
    g_obs,
    mono_params,
    dilution_rate,
    s_in,
    g_in,
):
    a_at = float(params_array[0])
    a_oa = float(params_array[1])

    at_pred, oa_pred, s_pred, g_pred = simulate_segment(
        t_start=t_start,
        y0=y0,
        times_eval=times_eval,
        a_at=a_at,
        a_oa=a_oa,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    eps = 1e-12

    at_res = np.log10(np.clip(at_pred, eps, None)) - np.log10(
        np.clip(at_obs, eps, None)
    )
    oa_res = np.log10(np.clip(oa_pred, eps, None)) - np.log10(
        np.clip(oa_obs, eps, None)
    )
    s_res = np.log10(np.clip(s_pred, eps, None)) - np.log10(np.clip(s_obs, eps, None))
    g_res = np.log10(np.clip(g_pred, eps, None)) - np.log10(np.clip(g_obs, eps, None))

    return np.concatenate([at_res, oa_res, s_res, g_res])


def fit_segment(
    t_start,
    y0,
    times_eval,
    at_obs,
    oa_obs,
    s_obs,
    g_obs,
    mono_params,
    dilution_rate,
    s_in,
    g_in,
):
    result = least_squares(
        fun=residuals_segment,
        x0=np.array([0.5, 0.5]),
        bounds=(np.array([1e-6, 1e-6]), np.array([1 - 1e-6, 1 - 1e-6])),
        args=(
            t_start,
            y0,
            times_eval,
            at_obs,
            oa_obs,
            s_obs,
            g_obs,
            mono_params,
            dilution_rate,
            s_in,
            g_in,
        ),
    )

    a_at = float(result.x[0])
    a_oa = float(result.x[1])

    at_pred, oa_pred, s_pred, g_pred = simulate_segment(
        t_start=t_start,
        y0=y0,
        times_eval=times_eval,
        a_at=a_at,
        a_oa=a_oa,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    return {
        "a_at": a_at,
        "a_oa": a_oa,
        "success": bool(result.success),
        "message": result.message,
        "cost": float(result.cost),
        "at_pred": at_pred,
        "oa_pred": oa_pred,
        "s_pred": s_pred,
        "g_pred": g_pred,
    }


def fit_single_phase_allocation(coculture_data, mono_params, dilution_rate):
    times = coculture_data["times"]
    at_obs = coculture_data["at_obs"]
    oa_obs = coculture_data["oa_obs"]
    s_obs = coculture_data["succinate_obs"]
    g_obs = coculture_data["glucose_obs"]

    y0 = [
        float(at_obs[0]),
        float(oa_obs[0]),
        float(s_obs[0]),
        float(g_obs[0]),
    ]
    s_in = float(s_obs[0])
    g_in = float(g_obs[0])

    single_fit = fit_segment(
        t_start=times[0],
        y0=y0,
        times_eval=times,
        at_obs=at_obs,
        oa_obs=oa_obs,
        s_obs=s_obs,
        g_obs=g_obs,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    t_sim = np.linspace(times[0], times[-1], 400)
    at_sim, oa_sim, s_sim, g_sim = simulate_segment(
        t_start=times[0],
        y0=y0,
        times_eval=t_sim,
        a_at=single_fit["a_at"],
        a_oa=single_fit["a_oa"],
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    return {
        "mode": "single",
        "times": times,
        "at_obs": at_obs,
        "oa_obs": oa_obs,
        "succinate_obs": s_obs,
        "glucose_obs": g_obs,
        "a_at": single_fit["a_at"],
        "a_oa": single_fit["a_oa"],
        "success": single_fit["success"],
        "message": single_fit["message"],
        "cost": single_fit["cost"],
        "s_in": s_in,
        "g_in": g_in,
        "t_sim": t_sim,
        "at_sim": at_sim,
        "oa_sim": oa_sim,
        "succinate_sim": s_sim,
        "glucose_sim": g_sim,
    }


def fit_two_phase_allocation(coculture_data, mono_params, dilution_rate, switch_time):
    times = coculture_data["times"]
    at_obs = coculture_data["at_obs"]
    oa_obs = coculture_data["oa_obs"]
    s_obs = coculture_data["succinate_obs"]
    g_obs = coculture_data["glucose_obs"]

    if switch_time <= times[0] or switch_time >= times[-1]:
        raise ValueError("SWITCH_TIME must lie inside the fitted time range.")

    pre_mask = times <= switch_time
    post_mask = times > switch_time

    if pre_mask.sum() < 2:
        raise ValueError(
            "Not enough pre-switch timepoints to fit pre-switch allocation."
        )
    if post_mask.sum() < 1:
        raise ValueError(
            "No post-switch timepoints available to fit post-switch allocation."
        )

    times_pre = times[pre_mask]
    times_post = times[post_mask]

    at_pre = at_obs[pre_mask]
    oa_pre = oa_obs[pre_mask]
    s_pre = s_obs[pre_mask]
    g_pre = g_obs[pre_mask]

    at_post = at_obs[post_mask]
    oa_post = oa_obs[post_mask]
    s_post = s_obs[post_mask]
    g_post = g_obs[post_mask]

    y0 = [
        float(at_obs[0]),
        float(oa_obs[0]),
        float(s_obs[0]),
        float(g_obs[0]),
    ]
    s_in = float(s_obs[0])
    g_in = float(g_obs[0])

    pre_fit = fit_segment(
        t_start=times_pre[0],
        y0=y0,
        times_eval=times_pre,
        at_obs=at_pre,
        oa_obs=oa_pre,
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
        a_at=pre_fit["a_at"],
        a_oa=pre_fit["a_oa"],
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    post_fit = fit_segment(
        t_start=switch_time,
        y0=state_switch,
        times_eval=times_post,
        at_obs=at_post,
        oa_obs=oa_post,
        s_obs=s_post,
        g_obs=g_post,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    t_pre_dense = np.linspace(times[0], switch_time, 250)
    t_post_dense = np.linspace(switch_time, times[-1], 250)

    at_pre_sim, oa_pre_sim, s_pre_sim, g_pre_sim = simulate_segment(
        t_start=times[0],
        y0=y0,
        times_eval=t_pre_dense,
        a_at=pre_fit["a_at"],
        a_oa=pre_fit["a_oa"],
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    state_switch_dense = [
        at_pre_sim[-1],
        oa_pre_sim[-1],
        s_pre_sim[-1],
        g_pre_sim[-1],
    ]

    at_post_sim, oa_post_sim, s_post_sim, g_post_sim = simulate_segment(
        t_start=switch_time,
        y0=state_switch_dense,
        times_eval=t_post_dense,
        a_at=post_fit["a_at"],
        a_oa=post_fit["a_oa"],
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    t_sim = np.concatenate([t_pre_dense, t_post_dense[1:]])
    at_sim = np.concatenate([at_pre_sim, at_post_sim[1:]])
    oa_sim = np.concatenate([oa_pre_sim, oa_post_sim[1:]])
    s_sim = np.concatenate([s_pre_sim, s_post_sim[1:]])
    g_sim = np.concatenate([g_pre_sim, g_post_sim[1:]])

    return {
        "mode": "switch",
        "times": times,
        "at_obs": at_obs,
        "oa_obs": oa_obs,
        "succinate_obs": s_obs,
        "glucose_obs": g_obs,
        "a_at_pre": pre_fit["a_at"],
        "a_oa_pre": pre_fit["a_oa"],
        "a_at_post": post_fit["a_at"],
        "a_oa_post": post_fit["a_oa"],
        "pre_success": pre_fit["success"],
        "post_success": post_fit["success"],
        "pre_message": pre_fit["message"],
        "post_message": post_fit["message"],
        "pre_cost": pre_fit["cost"],
        "post_cost": post_fit["cost"],
        "switch_time": switch_time,
        "s_in": s_in,
        "g_in": g_in,
        "t_sim": t_sim,
        "at_sim": at_sim,
        "oa_sim": oa_sim,
        "succinate_sim": s_sim,
        "glucose_sim": g_sim,
    }


def simulate_manual_allocation(coculture_data, mono_params, dilution_rate, a_at, a_oa):
    times = coculture_data["times"]
    at_obs = coculture_data["at_obs"]
    oa_obs = coculture_data["oa_obs"]
    s_obs = coculture_data["succinate_obs"]
    g_obs = coculture_data["glucose_obs"]

    y0 = [
        float(at_obs[0]),
        float(oa_obs[0]),
        float(s_obs[0]),
        float(g_obs[0]),
    ]
    s_in = float(s_obs[0])
    g_in = float(g_obs[0])

    t_sim = np.linspace(times[0], times[-1], 400)
    at_sim, oa_sim, s_sim, g_sim = simulate_segment(
        t_start=times[0],
        y0=y0,
        times_eval=t_sim,
        a_at=a_at,
        a_oa=a_oa,
        mono_params=mono_params,
        dilution_rate=dilution_rate,
        s_in=s_in,
        g_in=g_in,
    )

    return {
        "mode": "manual",
        "times": times,
        "at_obs": at_obs,
        "oa_obs": oa_obs,
        "succinate_obs": s_obs,
        "glucose_obs": g_obs,
        "a_at": a_at,
        "a_oa": a_oa,
        "s_in": s_in,
        "g_in": g_in,
        "t_sim": t_sim,
        "at_sim": at_sim,
        "oa_sim": oa_sim,
        "succinate_sim": s_sim,
        "glucose_sim": g_sim,
    }


def make_cfu_plot(fit):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["at_obs"],
            mode="markers",
            name="At data",
            marker=dict(
                color=colors["At"],
                size=8,
                opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=fit["t_sim"],
            y=fit["at_sim"],
            mode="lines",
            name="At fit",
            line=dict(color=colors["At"], width=3),
            opacity=0.8,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["oa_obs"],
            mode="markers",
            name="Oa data",
            marker=dict(
                color=colors["Oa"],
                size=8,
                opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=fit["t_sim"],
            y=fit["oa_sim"],
            mode="lines",
            name="Oa fit",
            line=dict(color=colors["Oa"], width=3),
            opacity=0.8,
        )
    )

    if fit["mode"] == "switch":
        fig.add_vline(
            x=fit["switch_time"],
            line=dict(color="gray", dash="dash", width=1.5),
        )

    fig.update_layout(
        width=200,
        height=180,
        title="At Oa Succinate Glucose",
        showlegend=False,
        xaxis=dict(title="Time (h)", range=[0, TIME_MAX], dtick=20),
        yaxis=dict(
            title="CFU/mL",
            type="log",
            exponentformat="power",
            range=[4, 10],
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

    if fit["mode"] == "switch":
        fig.add_vline(
            x=fit["switch_time"],
            line=dict(color="gray", dash="dash", width=1.5),
        )

    fig.update_layout(
        width=260,
        height=200,
        title="At Oa Succinate Glucose",
        showlegend=False,
        xaxis=dict(title="Time (h)", range=[0, TIME_MAX], dtick=20),
        yaxis=dict(
            title="Area",
            type="log",
            exponentformat="power",
        ),
    )

    fig = style_plot(fig, marker_size=8, line_thickness=3, font_size=11)
    return fig


def flatten_mono_params(mono_params):
    return {
        "mumax_at_succinate": mono_params["at"]["succinate"]["mumax"],
        "km_at_succinate": mono_params["at"]["succinate"]["km"],
        "yield_at_succinate": mono_params["at"]["succinate"]["yield"],
        "mumax_at_glucose": mono_params["at"]["glucose"]["mumax"],
        "km_at_glucose": mono_params["at"]["glucose"]["km"],
        "yield_at_glucose": mono_params["at"]["glucose"]["yield"],
        "mumax_oa_succinate": mono_params["oa"]["succinate"]["mumax"],
        "km_oa_succinate": mono_params["oa"]["succinate"]["km"],
        "yield_oa_succinate": mono_params["oa"]["succinate"]["yield"],
        "mumax_oa_glucose": mono_params["oa"]["glucose"]["mumax"],
        "km_oa_glucose": mono_params["oa"]["glucose"]["km"],
        "yield_oa_glucose": mono_params["oa"]["glucose"]["yield"],
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    output_paths = get_output_paths(USE_SWITCH)

    mono_params = load_mono_parameters(MONO_PARAMETER_CSV)
    resource_df = pd.read_csv("normalized_data.csv", index_col="id")
    cfu_df = pd.read_csv("cfus.csv")

    coculture_data = prepare_coculture_data(resource_df=resource_df, cfu_df=cfu_df)

    if MANUAL_A:
        fit = simulate_manual_allocation(
            coculture_data=coculture_data,
            mono_params=mono_params,
            dilution_rate=DILUTION_RATE,
            a_at=A_AT,
            a_oa=A_OA,
        )

        result_row = {
            "mode": "manual",
            "time_max": TIME_MAX,
            "dilution_rate": DILUTION_RATE,
            "a_at_succinate": fit["a_at"],
            "a_at_glucose": 1 - fit["a_at"],
            "a_oa_succinate": fit["a_oa"],
            "a_oa_glucose": 1 - fit["a_oa"],
            **flatten_mono_params(mono_params),
        }

        print(f"\nManual At succinate allocation = {fit['a_at']:.6f}")
        print(f"Manual At glucose allocation   = {1 - fit['a_at']:.6f}")
        print(f"Manual Oa succinate allocation = {fit['a_oa']:.6f}")
        print(f"Manual Oa glucose allocation   = {1 - fit['a_oa']:.6f}")

    elif USE_SWITCH:
        fit = fit_two_phase_allocation(
            coculture_data=coculture_data,
            mono_params=mono_params,
            dilution_rate=DILUTION_RATE,
            switch_time=SWITCH_TIME,
        )

        result_row = {
            "mode": "switch",
            "switch_time": fit["switch_time"],
            "time_max": TIME_MAX,
            "dilution_rate": DILUTION_RATE,
            "a_at_pre_succinate": fit["a_at_pre"],
            "a_at_pre_glucose": 1 - fit["a_at_pre"],
            "a_oa_pre_succinate": fit["a_oa_pre"],
            "a_oa_pre_glucose": 1 - fit["a_oa_pre"],
            "a_at_post_succinate": fit["a_at_post"],
            "a_at_post_glucose": 1 - fit["a_at_post"],
            "a_oa_post_succinate": fit["a_oa_post"],
            "a_oa_post_glucose": 1 - fit["a_oa_post"],
            "pre_success": fit["pre_success"],
            "post_success": fit["post_success"],
            "pre_cost": fit["pre_cost"],
            "post_cost": fit["post_cost"],
            "pre_message": fit["pre_message"],
            "post_message": fit["post_message"],
            **flatten_mono_params(mono_params),
        }

    else:
        fit = fit_single_phase_allocation(
            coculture_data=coculture_data,
            mono_params=mono_params,
            dilution_rate=DILUTION_RATE,
        )

        result_row = {
            "mode": "single",
            "time_max": TIME_MAX,
            "dilution_rate": DILUTION_RATE,
            "a_at_succinate": fit["a_at"],
            "a_at_glucose": 1 - fit["a_at"],
            "a_oa_succinate": fit["a_oa"],
            "a_oa_glucose": 1 - fit["a_oa"],
            "success": fit["success"],
            "cost": fit["cost"],
            "message": fit["message"],
            **flatten_mono_params(mono_params),
        }

    results_df = pd.DataFrame([result_row])
    results_df.to_csv(output_paths["result_csv"], index=False)

    fig_cfu = make_cfu_plot(fit)
    fig_cfu.write_image(os.path.join(OUTPUT_DIR, output_paths["cfu_fname"]))

    fig_resource = make_resource_plot(fit)
    fig_resource.write_image(os.path.join(OUTPUT_DIR, output_paths["resource_fname"]))

    print("\n" + results_df.to_string(index=False))


if __name__ == "__main__":
    main()
