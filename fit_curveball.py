import curveball
import pandas as pd
from style import *
import plotly.graph_objects as go
from experiment import Species, Experiment
import numpy as np
from plotly.subplots import make_subplots

e = Experiment(d="data/251018_succinate_glucose_plate_reader/metaod/")
carbon_sources = ["Succinate", "Glucose"]
species = ["Oa", "At"]

figures = {
    "Oa": make_subplots(rows=1, cols=2, subplot_titles=carbon_sources),
    "At": make_subplots(rows=1, cols=2, subplot_titles=carbon_sources),
}
models = []


def filter_time(xs, ys, x0=0, x1=27):
    xs_f = []
    ys_f = []
    for xi, yi in zip(xs, ys):
        mask = (xi >= x0) & (xi <= x1)
        xs = xi[mask]
        ys = yi[mask]
        xs_f.append(xs)
        ys_f.append(ys)
    return xs_f, ys_f


for sp in species:
    for cs in carbon_sources:
        c = e.get_condition(
            [sp],
            cs,
            15,
            "OD",
        )
        dfs = []
        xs, ys = c.xs, c.ys
        xs, ys = filter_time(xs, ys, x0=0, x1=27)
        for i in range(3):
            df = pd.DataFrame(
                {"Time": xs[i], "OD": ys[i], "Well": str(i), "Strain": sp}
            )
            dfs.append(df)

        m = curveball.models.fit_model(
            pd.concat(dfs),
            PLOT=False,
            PRINT=False,
            param_guess={"y0": np.average([y[0] for y in ys])},
            param_fix=["y0"],
        )
        m = sorted(m, key=lambda r: r.model.name)
        models.append(m)

for i, m in enumerate(models[:2]):
    fig = figures[species[0]]
    cs = carbon_sources[i]
    fit = m[2]
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
                line=dict(
                    width=1.2,
                ),
            ),
        ),
        row=1,
        col=i + 1,
    )
    fig.add_trace(
        go.Scatter(
            x=fit.userkws["t"],
            y=fit.best_fit,
            mode="lines",
            name="Fit",
            line=dict(
                color=colors["Succinate+Glucose Outflow"],
            ),
        ),
        row=1,
        col=i + 1,
    )
    fig.add_annotation(
        xref="x domain",
        yref="y domain",
        x=0.01,
        y=0.99,
        text=f"Oa μmax: {fit.params['r'].value:.3f} h⁻¹",
        row=1,
        col=i + 1,
        showarrow=False,
    )
    fig.update_layout(showlegend=False)

for i, m in enumerate(models[2:4]):
    fig = figures[species[1]]
    cs = carbon_sources[i]
    fit = m[2]
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
                line=dict(
                    width=1.2,
                ),
            ),
        ),
        row=1,
        col=i + 1,
    )
    fig.add_trace(
        go.Scatter(
            x=fit.userkws["t"],
            y=fit.best_fit,
            mode="lines",
            name="Fit",
            line=dict(
                color=colors["Succinate+Glucose Outflow"],
            ),
        ),
        row=1,
        col=i + 1,
    )
    fig.add_annotation(
        xref="x domain",
        yref="y domain",
        x=0.01,
        y=0.99,
        text=f"At μmax: {fit.params['r'].value:.3f} h⁻¹",
        row=1,
        col=i + 1,
        showarrow=False,
    )
    fig.update_layout(showlegend=False)


fnames = ["plots/fitting/oa_curveball.svg", "plots/fitting/at_curveball.svg"]
for fig, fname in zip(figures.values(), fnames):
    fig = style_plot(fig, marker_size=5, line_thickness=3, font_size=12)
    fig.write_image(fname)
