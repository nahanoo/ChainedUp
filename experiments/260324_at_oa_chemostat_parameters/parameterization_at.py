import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import conditions as cond
from style import *


DILUTION_RATE = 0.3
OUTPUT_DIR = "plots/parameterization"
PARAMETER_CSV = "parameterization_results.csv"

AT_GLUCOSE_MUMAX = 0.696  # set this to the fitted mu_max for At on glucose


SPECS = {
    "oa_succinate": {
        "species": "oa",
        "wells": cond.OA_SUCCINATE,
        "times": np.array(cond.OA_TIMEPOINTS, dtype=float),
        "substrate_col": "succinate_area_normalized",
        "reactor": "Succinate",
        "mumax": 0.877,
        "label": "Succinate",
        "color_key": "Succinate",
        "cfu_title": "Oa Succinate",
        "resource_title": "Oa Succinate",
        "cfu_fname": "oa_succinate_cfu_fit.svg",
        "resource_fname": "oa_succinate_resource_fit.svg",
    },
    "oa_glucose": {
        "species": "oa",
        "wells": cond.OA_GLUCOSE,
        "times": np.array(cond.OA_TIMEPOINTS, dtype=float),
        "substrate_col": "glucose_area_normalized",
        "reactor": "Glucose",
        "mumax": 0.573,
        "label": "Glucose",
        "color_key": "Glucose",
        "cfu_title": "Oa Glucose",
        "resource_title": "Oa Glucose",
        "cfu_fname": "oa_glucose_cfu_fit.svg",
        "resource_fname": "oa_glucose_resource_fit.svg",
    },
    "at_glucose": {
        "species": "at",
        "wells": cond.AT_SUCCINATE_GLUCOSE,
        "times": np.array(cond.AT_TIMEPOINTS, dtype=float),
        "substrate_col": "glucose_area_normalized",
        "reactor": "Glucose",
        "mumax": AT_GLUCOSE_MUMAX,
        "label": "Glucose",
        "color_key": "Glucose",
        "cfu_title": "At Glucose",
        "resource_title": "At Glucose",
        "cfu_fname": "at_glucose_cfu_fit.svg",
        "resource_fname": "at_glucose_resource_fit.svg",
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


def prepare_condition_data(
    resource_df,
    cfu_df,
    wells,
    times,
    substrate_col,
    reactor,
    species,
):
    substrate_df = pd.DataFrame(
        {
            "time": times,
            "s_obs": [resource_df.loc[well, substrate_col] for well in wells],
        }
    )

    biomass_df = (
        cfu_df[
            (cfu_df["species"].str.lower() == species.lower())
            & (cfu_df["reactor"] == reactor)
        ]
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
        raise ValueError(
            f"No overlapping timepoints found for species={species!r}, reactor={reactor!r}."
        )

    return (
        merged["time"].to_numpy(dtype=float),
        merged["n_obs"].to_numpy(dtype=float),
        merged["s_obs"].to_numpy(dtype=float),
    )


def residuals(log_params, times, n_obs, s_obs, mumax, dilution_rate, s_in):
    km, yld = 10**log_params

    n0 = n_obs[0]
    s0 = s_obs[0]

    n_pred, s_pred = simulate(
        times=times,
        n0=n0,
        s0=s0,
        mumax=mumax,
        km=km,
        yld=yld,
        dilution_rate=dilution_rate,
        s_in=s_in,
    )

    eps = 1e-12

    n_res = np.log10(n_pred + eps) - np.log10(n_obs + eps)
    s_res = np.log10(s_pred + eps) - np.log10(s_obs + eps)

    return np.concatenate([n_res, s_res])


def fit_condition(times, n_obs, s_obs, mumax, dilution_rate):
    s_in = s_obs[0]

    x0 = np.log10([1e-2, 1e8])
    lower = np.log10([1e-8, 1e4])
    upper = np.log10([1e2, 1e12])

    result = least_squares(
        fun=residuals,
        x0=x0,
        bounds=(lower, upper),
        args=(times, n_obs, s_obs, mumax, dilution_rate, s_in),
    )

    km, yld = 10**result.x

    n_pred, s_pred = simulate(
        times=times,
        n0=n_obs[0],
        s0=s_obs[0],
        mumax=mumax,
        km=km,
        yld=yld,
        dilution_rate=dilution_rate,
        s_in=s_in,
    )

    return {
        "km": km,
        "yield": yld,
        "success": result.success,
        "message": result.message,
        "cost": result.cost,
        "times": times,
        "n_obs": n_obs,
        "s_obs": s_obs,
        "n_pred": n_pred,
        "s_pred": s_pred,
        "n0": n_obs[0],
        "s0": s_obs[0],
        "s_in": s_in,
    }


def make_cfu_fit_plot(fit, spec):
    color = colors[spec["color_key"]]

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["n_obs"],
            mode="markers",
            name="Data",
            marker=dict(
                color=color,
                size=8,
                opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["n_pred"],
            mode="lines",
            name="Fit",
            line=dict(color=color, width=3),
            opacity=0.8,
        )
    )

    fig.update_layout(
        width=200,
        height=180,
        title=spec["cfu_title"],
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95,
            bgcolor="rgba(255,255,255,0)",
        ),
        xaxis=dict(title="Time (h)"),
        yaxis=dict(
            title="CFU/mL",
            type="log",
            exponentformat="power",
            range=[6, 10],
        ),
    )

    fig = style_plot(fig, marker_size=8, line_thickness=3, font_size=11)
    return fig


def make_resource_fit_plot(fit, spec):
    color = colors[spec["color_key"]]

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["s_obs"],
            mode="markers",
            name="Data",
            marker=dict(
                color=color,
                size=8,
                opacity=0.8,
                line=dict(color="black", width=1.2),
            ),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=fit["times"],
            y=fit["s_pred"],
            mode="lines",
            name="Fit",
            line=dict(color=color, width=3),
            opacity=0.8,
        )
    )

    fig.update_layout(
        width=200,
        height=180,
        title=spec["resource_title"],
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95,
            bgcolor="rgba(255,255,255,0)",
        ),
        xaxis=dict(title="Time (h)"),
        yaxis=dict(
            title="Area",
            type="log",
            exponentformat="power",
        ),
    )

    fig = style_plot(fig, marker_size=8, line_thickness=3, font_size=11)
    return fig


def validate_specs(specs):
    missing_mumax = [name for name, spec in specs.items() if spec["mumax"] is None]
    if missing_mumax:
        raise ValueError(
            "Missing mumax for the following conditions: "
            + ", ".join(missing_mumax)
            + ". Set AT_GLUCOSE_MUMAX before running."
        )


def main():
    validate_specs(SPECS)
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
            species=spec["species"],
        )

        fit = fit_condition(
            times=times,
            n_obs=n_obs,
            s_obs=s_obs,
            mumax=spec["mumax"],
            dilution_rate=DILUTION_RATE,
        )

        results.append(
            {
                "condition": name,
                "species": spec["species"],
                "reactor": spec["reactor"],
                "substrate": spec["label"],
                "mumax": spec["mumax"],
                "Km": fit["km"],
                "yield": fit["yield"],
                "n0": fit["n0"],
                "s0": fit["s0"],
                "s_in": fit["s_in"],
                "cost": fit["cost"],
                "success": fit["success"],
                "message": fit["message"],
            }
        )

        fig_cfu = make_cfu_fit_plot(fit=fit, spec=spec)
        fig_cfu.write_image(os.path.join(OUTPUT_DIR, spec["cfu_fname"]))

        fig_resource = make_resource_fit_plot(fit=fit, spec=spec)
        fig_resource.write_image(os.path.join(OUTPUT_DIR, spec["resource_fname"]))

    results_df = pd.DataFrame(results)
    results_df.to_csv(PARAMETER_CSV, index=False)

    print(results_df.to_string(index=False))


if __name__ == "__main__":
    main()
