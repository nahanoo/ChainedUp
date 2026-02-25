"""
Consumer–resource model with Baranyi lag and Monod uptake.

Workflow:
1) Fit growth curve with curveball (logistic-lag) to get r and q0 (and Kcarry).
2) Convert Kcarry + resource concentration to a yield Y.
3) Fit ONLY Km (Monod half-saturation constant) by matching the ODE model to OD(t).

Units:
- t in hours
- resource R in mM
- Km in mM (bounds are set in µM->mM)
- N is OD (or biomass proxy)
- Y is OD per mM resource
"""

import os.path as path

import lmfit
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.integrate import solve_ivp

from experiment import Experiment
from fit_curveball import condition_to_df, fit_curveball
from style import *


def sort_unique_xy(t, y, how="mean"):
    """Sort (t,y) by time and deduplicate t (needed for solve_ivp t_eval)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(t) & np.isfinite(y)
    t, y = t[m], y[m]

    order = np.argsort(t)
    t, y = t[order], y[order]

    if how == "first":
        t_u, idx = np.unique(t, return_index=True)
        return t_u, y[idx]

    if how == "mean":
        t_u = np.unique(t)
        y_u = np.array([y[t == tu].mean() for tu in t_u], dtype=float)
        return t_u, y_u

    raise ValueError("how must be 'first' or 'mean'")


def alpha_baranyi(t, q0, v):
    """
    Baranyi adjustment function alpha(t) in [0,1].

    Convention here:
    - q0 <= 0 (or None) disables lag (alpha = 1).
    """
    if q0 is None or q0 <= 0:
        return 1.0
    return q0 / (q0 + np.exp(-v * t))


def cr_lag_rhs(t, y, r, Km, q0, Y, v, D=0.0, M=0.0):
    """ODE RHS for consumer (N) and resource (R)."""
    N, R = y
    R = max(R, 0.0)
    N = max(N, 0.0)

    a = alpha_baranyi(t, q0=q0, v=v)
    mu = r * R / (R + Km)

    dNdt = a * mu * N - D * N
    dRdt = -(1.0 / Y) * a * mu * N - D * R + D * M
    return [dNdt, dRdt]


def simulate_cr_lag(t, N0, R0, r, Km, q0, Y, v=None, method="LSODA", D=0.0, M=0.0):
    """
    Simulate model at times t. If v is not provided, we set v=r (common constraint).
    """
    if v is None:
        v = r

    t = np.asarray(t, dtype=float)
    if t.size < 2:
        raise ValueError("Need at least 2 timepoints.")

    # solve_ivp requires strictly increasing t_eval
    if not np.all(np.diff(t) > 0):
        raise ValueError("t must be strictly increasing (use sort_unique_xy).")

    N0 = float(max(N0, 0.0))
    R0 = float(max(R0, 0.0))

    sol = solve_ivp(
        fun=lambda tt, yy: cr_lag_rhs(tt, yy, r=r, Km=Km, q0=q0, Y=Y, v=v, D=D, M=M),
        t_span=(float(t[0]), float(t[-1])),
        y0=[N0, R0],
        t_eval=t,
        method=method,
        rtol=1e-7,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    N_hat, R_hat = sol.y
    return N_hat, R_hat


def fit_Km_lmfit(
    t,
    N_data,
    *,
    N0,
    R0,
    r,
    q0,
    Y,
    v=None,
    Km_guess=0.1,  # mM
    Km_min_uM=0.001,  # µM
    Km_max_mM=20.0,  # mM
    dedup="mean",
    method="least_squares",
    D=0.0,
    M=0.0,
):
    """
    Fit ONLY Km (Monod half-saturation constant).
    Bounds:
      Km ∈ [Km_min_uM, Km_max_mM] converted to mM for the model.
    """
    t, N_data = sort_unique_xy(t, N_data, how=dedup)
    Km_min_mM = Km_min_uM * 1e-3  # µM -> mM

    params = lmfit.Parameters()
    params.add("Km", value=Km_guess, min=Km_min_mM, max=Km_max_mM)

    def residual(p):
        Km = p["Km"].value
        N_hat, _ = simulate_cr_lag(
            t,
            N0=N0,
            R0=R0,
            r=r,
            Km=Km,
            q0=q0,
            Y=Y,
            v=v,
            D=D,
            M=M,
        )
        return N_hat - N_data

    return lmfit.minimize(residual, params, method=method)


def rolling_log_slope(t, N, *, window_points=7):
    """
    Rolling slope of log(N) vs t using a sliding linear regression window.
    Returns (t_center, slope) where slope is d/dt log(N) [h^-1].

    Assumes t is strictly increasing and N > 0.
    """
    t = np.asarray(t, dtype=float)
    N = np.asarray(N, dtype=float)

    if t.size < window_points:
        raise ValueError("Not enough points for rolling slope window.")

    lnN = np.log(N)

    ts = []
    slopes = []
    for i in range(window_points, t.size + 1):
        tw = t[i - window_points : i]
        yw = lnN[i - window_points : i]

        A = np.vstack([tw, np.ones_like(tw)]).T
        b, a0 = np.linalg.lstsq(A, yw, rcond=None)[0]  # slope, intercept

        ts.append(tw.mean())
        slopes.append(b)

    return np.asarray(ts), np.asarray(slopes)


def estimate_rmax_chemostat(
    t,
    N,
    *,
    D,
    q0=0.0,
    v=None,
    window_points=7,
    min_od=1e-4,
    dedup="mean",
):
    """
    Estimate an *effective* maximum specific growth rate from chemostat OD(t).

    Uses:
      d/dt log(N) = alpha(t) * mu(t) - D
      => mu_hat(t) = (d/dt log(N) + D) / alpha(t)

    Then returns:
      r_hat = max_t mu_hat(t)

    Notes:
    - mu(t) <= r (Monod), so r_hat is typically a LOWER BOUND on r unless early-time
      resource is saturating (R >> Km).
    - If q0 <= 0, alpha(t)=1 (no lag correction).
    """
    t, N = sort_unique_xy(t, N, how=dedup)

    N = np.asarray(N, dtype=float)
    t = np.asarray(t, dtype=float)

    # baseline shift if OD has negatives (common with blank subtraction)
    nmin = np.nanmin(N)
    if not np.isfinite(nmin):
        raise ValueError("N contains no finite values.")
    if nmin <= 0:
        N = N - nmin + min_od

    m = np.isfinite(N) & (N > min_od)
    t2 = t[m]
    N2 = N[m]

    if t2.size < max(window_points, 3):
        raise ValueError(
            "Not enough positive OD points to estimate r from chemostat data."
        )

    if not np.all(np.diff(t2) > 0):
        raise ValueError(
            "t must be strictly increasing after filtering (check duplicates)."
        )

    if v is None:
        # only used when q0>0; for a first-pass estimate this is fine
        v = 1.0

    tt, slopes = rolling_log_slope(t2, N2, window_points=window_points)

    alpha = np.array([alpha_baranyi(ti, q0=q0, v=v) for ti in tt], dtype=float)
    alpha = np.clip(alpha, 1e-9, 1.0)

    mu_hat = (slopes + float(D)) / alpha

    i = int(np.nanargmax(mu_hat))
    r_hat = float(mu_hat[i])

    info = dict(
        t=tt,
        slope_net=slopes,
        mu_hat=mu_hat,
        t_max=float(tt[i]),
        mu_max=float(mu_hat[i]),
    )
    return r_hat, info


def plot_fit(t_raw, y_raw, t_fit, y_fit, r, Km, fig, *, label=None):
    fig.add_trace(
        go.Scatter(
            x=t_raw,
            y=y_raw,
            mode="markers",
            name="Data" if label is None else f"Data ({label})",
            marker=dict(
                symbol="circle",
                color=colors["Glucose"],
                opacity=0.55,
                line=dict(width=1.2),
            ),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t_fit,
            y=y_fit,
            mode="lines",
            name="Model" if label is None else f"Model ({label})",
            line=dict(color=colors["Succinate+Glucose Outflow"]),
            showlegend=False,
        )
    )
    fig.add_annotation(
        x=t_fit[-1],
        y=y_fit[-1],
        text=f"r {r:.3f} h⁻¹, Km {Km:.4f} mM",
        showarrow=False,
        yshift=10,
    )
    return fig


def fit_Km_chemostat(
    csv_path,
    *,
    species="Ct",
    substrate="acetate",
    replicate="replicate_1",
    r=0.4,
    q0=0.0,  # q0<=0 disables lag in alpha_baranyi()
    Y=0.03,
    D=0.1,
    M=7.5,  # inflow resource concentration [mM]
    Km_guess=0.001,
    R0=None,  # if None -> start at inflow M
    r_from_data=False,
):
    csv_path = path.expanduser(csv_path)
    df = pd.read_csv(csv_path)

    df = df[
        (df["species"] == species)
        & (df["substrate"] == substrate)
        & (df["replicate"] == replicate)
    ].copy()

    t_raw = df["time"].to_numpy()
    y_raw = df["OD"].to_numpy()

    if r_from_data:
        r_est, info = estimate_rmax_chemostat(
            t_raw, y_raw, D=D, q0=q0, window_points=60 * 5
        )
        r = r_est

    # ensure solve_ivp is happy
    t_fit, y_fit_data = sort_unique_xy(t_raw, y_raw, how="mean")

    N0 = float(max(y_fit_data[0], 0.08))
    if R0 is None:
        R0 = M

    res = fit_Km_lmfit(
        t_fit,
        y_fit_data,
        N0=N0,
        R0=R0,
        r=r,
        q0=q0,
        Y=Y,
        Km_guess=Km_guess,
        D=D,
        M=M,
    )

    Km = res.params["Km"].value
    y_hat = simulate_cr_lag(t_fit, N0=N0, R0=R0, r=r, Km=Km, q0=q0, Y=Y, D=D, M=M)[0]

    fig = go.Figure()
    plot_fit(t_raw, y_raw, t_fit, y_hat, r, Km, fig)
    fig = style_plot(fig, marker_size=6, line_thickness=3, font_size=12)

    if r_from_data:
        fig.add_annotation(
            x=t_fit[0],
            y=np.nanmax(y_raw),
            text=f"r_est from chemostat: {r:.3f} h⁻¹",
            showarrow=False,
            xanchor="left",
            yanchor="top",
        )

    return res, fig


def fit_Km_batch():
    species = [["At"], ["Oa"], ["At", "Oa"]]
    carbon_sources = ["Succinate", "Glucose"]
    e = Experiment(d="data/251018_succinate_glucose_plate_reader/metaod/")
    concs = [5, 7.5, 15, 30]

    for sp in species:
        for carbon_source in carbon_sources:
            fig = go.Figure()

            for conc in concs:
                cond = e.get_condition(sp, carbon_source, conc, "OD")
                df, ys = condition_to_df(cond, sp, x0=0, x1=27)
                models = fit_curveball(df, ys)
                fit = models[2]  # logistic-lag model
                p = fit.params

                q0 = p["q0"].value
                r = p["r"].value
                N0 = p["y0"].value
                Y = (p["K"].value - p["y0"].value) / conc

                t_raw, y_raw = fit.userkws["t"], fit.data

                res = fit_Km_lmfit(
                    t_raw,
                    y_raw,
                    N0=N0,
                    R0=conc,
                    r=r,
                    q0=q0,
                    Y=Y,
                    Km_guess=0.05,
                )
                lmfit.report_fit(res)

                t_fit, y_fit_data = sort_unique_xy(t_raw, y_raw, how="mean")
                Km = res.params["Km"].value
                y_hat = simulate_cr_lag(t_fit, N0=N0, R0=conc, r=r, Km=Km, q0=q0, Y=Y)[
                    0
                ]

                plot_fit(t_raw, y_raw, t_fit, y_hat, r, Km, fig, label=f"{conc} mM")

            fig = style_plot(fig, marker_size=5, line_thickness=3, font_size=12)
            fname = f"plots/fitting/{'+'.join(sp)}_{carbon_source}.svg"
            fig.write_image(fname)


if __name__ == "__main__":
    # --- Chemostat example ---
    csv = "~/apps/curveball/data/chemostat_data/data.csv"

    # quick estimate only (no fitting)
    df0 = pd.read_csv(path.expanduser(csv))
    df0 = df0[
        (df0["species"] == "Oa")
        & (df0["substrate"] == "succinate")
        & (df0["replicate"] == "replicate_0")
    ]
    r_hat, info = estimate_rmax_chemostat(
        df0["time"].to_numpy(),
        df0["OD"].to_numpy(),
        D=0.3,
        q0=0,
        window_points=30,
        min_od=1e-4,
    )
    print(
        f"r_est (chemostat, from OD slopes): {r_hat:.4f} h^-1 at t≈{info['t_max']:.2f} h"
    )

    # then fit Km using either your provided r OR the estimated r
    res, fig = fit_Km_chemostat(
        csv,
        species="Oa",
        substrate="succinate",
        replicate="replicate_0",
        r=r_hat,
        q0=0,
        Y=0.066,
        D=0.3,
        M=7.5,
        Km_guess=0.001,
        r_from_data=True,  # <-- uses estimate_rmax_chemostat() internally
    )
    lmfit.report_fit(res)
    fig.write_image("tmp.svg")

    # --- Batch fitting ---
    # fit_Km_batch()
