from fit_curveball import condition_to_df, fit_curveball
from style import *
import numpy as np
from scipy.integrate import solve_ivp
import lmfit
from experiment import Experiment
import plotly.graph_objects as go


def sort_unique_xy(t, y, how="mean"):
    """
    Make t strictly increasing for solve_ivp(t_eval=...).
    Also keeps y aligned.

    how:
      - "first": keep first y at duplicated t
      - "mean": average y at duplicated t
    """
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
    return q0 / (q0 + np.exp(-v * t))


def cr_lag_rhs(t, y, r, Km, q0, Y, v):
    N, R = y

    # guard against tiny negatives from solver noise
    R = max(R, 0.0)
    N = max(N, 0.0)

    a = alpha_baranyi(t, q0=q0, v=v)
    mu = r * R / (R + Km)

    dNdt = a * mu * N
    dRdt = -(1.0 / Y) * a * mu * N
    return [dNdt, dRdt]


def simulate_cr_lag(t, N0, R0, r, Km, q0, Y, v=None, method="LSODA"):
    if v is None:
        v = r  # consistent with having only q0 from the logistic-lag style fit

    t = np.asarray(t, dtype=float)
    if t.size < 2:
        raise ValueError("Need at least 2 timepoints after sorting/dedup.")
    if not np.all(np.diff(t) > 0):
        raise ValueError("t must be strictly increasing (use sort_unique_xy).")

    sol = solve_ivp(
        fun=lambda tt, yy: cr_lag_rhs(tt, yy, r=r, Km=Km, q0=q0, Y=Y, v=v),
        t_span=(float(t[0]), float(t[-1])),
        y0=[float(N0), float(R0)],
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
):
    # solve_ivp needs sorted unique timepoints
    t, N_data = sort_unique_xy(t, N_data, how=dedup)

    Km_min_mM = Km_min_uM * 1e-3  # µM -> mM  (1 µM = 1e-3 mM)

    params = lmfit.Parameters()
    params.add("Km", value=Km_guess, min=Km_min_mM, max=Km_max_mM)

    def residual(p):
        Km = p["Km"].value
        N_hat, _ = simulate_cr_lag(t, N0=N0, R0=R0, r=r, Km=Km, q0=q0, Y=Y, v=v)
        return N_hat - N_data

    result = lmfit.minimize(residual, params, method=method)
    return result


# --- run ---
e = Experiment(d="data/251018_succinate_glucose_plate_reader/metaod/")
conc = 15
sp = "Oa"
cond = e.get_condition([sp], "Succinate", conc, "OD")
df, ys = condition_to_df(cond, sp, x0=0, x1=27)

models = fit_curveball(df, ys)
fit = models[2]
params = fit.params

q0, r, N0, Y = (
    params["q0"].value,
    params["r"].value,
    params["y0"].value,
    (params["K"].value - params["y0"].value) / conc,
    # your yield assumption: Y = Kcarry / R0
)

t = fit.userkws["t"]
od = fit.data

result = fit_Km_lmfit(t, od, N0=N0, R0=conc, r=r, q0=q0, Y=Y, Km_guess=0.05)
lmfit.report_fit(result)
t, N_data = sort_unique_xy(fit.userkws["t"], fit.data, how="mean")

y_fit = simulate_cr_lag(
    t, N0=N0, R0=conc, r=r, Km=result.params["Km"].value, q0=q0, Y=Y
)[0]
fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=fit.userkws["t"],
        y=fit.data,
        mode="markers",
        name="Data",
        marker=dict(
            symbol="circle",
            color=colors["Glucose"],
            opacity=0.55,
            line=dict(width=1.2),
        ),
    )
)
fig.add_trace(
    go.Scatter(
        x=t,
        y=y_fit,
        mode="lines",
        name="Fit",
        line=dict(color=colors["Succinate+Glucose Outflow"]),
    )
)
fig.add_annotation(
    xref="x domain",
    yref="y domain",
    x=0.01,
    y=0.99,
    text=f"{'At'} μmax: {fit.params['r'].value:.3f} h⁻¹<br>Km: {result.params['Km'].value:.3f} mM",
    showarrow=False,
)
fig = style_plot(fig, marker_size=5, line_thickness=3, font_size=12)
fig.write_image("tmp.svg")
