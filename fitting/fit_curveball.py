import curveball
import pandas as pd
from style import *
import plotly.graph_objects as go
from experiment import Experiment
import numpy as np
from plotly.subplots import make_subplots


def filter_time(xs, ys, x0=0, x1=27):
    xs_f, ys_f = [], []
    for xi, yi in zip(xs, ys):
        mask = (xi >= x0) & (xi <= x1)
        xs_f.append(xi[mask])
        ys_f.append(yi[mask])
    return xs_f, ys_f


def condition_to_df(cond, sp, x0=0, x1=27):
    xs, ys = filter_time(cond.xs, cond.ys, x0=x0, x1=x1)
    dfs = [
        pd.DataFrame(
            {"Time": xs[i], "OD": ys[i], "Well": str(i), "Strain": "+".join(sp)}
        )
        for i in range(3)
    ]
    return pd.concat(dfs), ys


def fit_curveball(df, ys):
    m = curveball.models.fit_model(
        df,
        PLOT=False,
        PRINT=False,
        param_guess={"y0": np.average([y[0] for y in ys])},
        param_fix=["y0"],
    )
    return sorted(m, key=lambda r: r.model.name)


def plot_fit(fig, fit, sp):
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
        ),
    )
    fig.add_trace(
        go.Scatter(
            x=fit.userkws["t"],
            y=fit.best_fit,
            mode="lines",
            name="Fit",
            line=dict(color=colors["Succinate+Glucose Outflow"]),
        ),
    )
    fig.add_annotation(
        xref="x domain",
        yref="y domain",
        x=0.01,
        y=0.99,
        text=f"μmax: {fit.params['r'].value:.3f} h⁻¹",
        showarrow=False,
    )
    fig.update_layout(showlegend=False)
    fig = style_plot(fig, marker_size=5, line_thickness=3, font_size=12)
    return fig


def fit_km(e, sp, carbon_sources, conc=15, metric="OD", x0=0, x1=27):
    for j, cs in enumerate(carbon_sources, start=1):
        cond = e.get_condition([sp], cs, conc, metric)
        df, ys = condition_to_df(cond, sp, x0=x0, x1=x1)
        models = fit_curveball(df, ys)
        fit = models[2]  # keep your original choice
    return fit
