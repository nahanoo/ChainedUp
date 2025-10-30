from dash import Dash, dcc, html, Output, Input
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from style import *
from plotly.subplots import make_subplots
from model import Model
from experiment import Species
from os import path

conc = 15
p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)
at = Species("At", params.loc["At"])
oa = Species("Oa", params.loc["Oa"])
global D
D = 0.1
M = Model(at, oa, None, np.linspace(0, 200, 500), conc, D)
# --- Dash app ---
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container(
    [
        # dbc.Row(html.H1("Fries Burger-Veggie model")),
        dbc.Row(html.H3("Parameters")),
        dbc.Row([dbc.Col([html.H4("At")], width=4), dbc.Col(html.H4("Oa"), width=4)]),
        # parameter controls
        dbc.Row(
            [
                dbc.Col(html.Span("μ Succinate"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="μ_at_succ",
                        type="number",
                        value=at.v_succ,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("μ Succinate"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="μ_oa_succ",
                        type="number",
                        value=oa.v_succ,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Span("μ Glucose"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="μ_at_gluc",
                        type="number",
                        value=at.v_gluc,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("μ Glucose"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="μ_oa_gluc",
                        type="number",
                        value=oa.v_gluc,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Span("K Succinate"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="at_K_succ",
                        type="number",
                        value=at.K_succ,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("K Succinate"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="oa_K_succ",
                        type="number",
                        value=oa.K_succ,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Span("K Glucose"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="at_K_gluc",
                        type="number",
                        value=at.K_gluc,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("K Glucose"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="oa_K_gluc",
                        type="number",
                        value=oa.K_gluc,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Span("Resource preference"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="at_a",
                        type="number",
                        value=at.a,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("Resource preference"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="oa_a",
                        type="number",
                        value=oa.a,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Span("Dilution rate"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="D",
                        type="number",
                        value=D,
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        # graphs
        dbc.Row(html.H4("At mono-culture")),
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id="at_succinate"), width=3),
                dbc.Col(dcc.Graph(id="at_glucose"), width=3),
                dbc.Col(dcc.Graph(id="at_succinate+glucose"), width=3),
            ]
        ),
        dbc.Row(html.H4("Oa mono-culture")),
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id="oa_succinate"), width=3),
                dbc.Col(dcc.Graph(id="oa_glucose"), width=3),
                dbc.Col(dcc.Graph(id="oa_succinate+glucose"), width=3),
            ]
        ),
        dbc.Row(html.H4("At + Oa co-culture")),
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id="at_oa_succinate"), width=3),
                dbc.Col(dcc.Graph(id="at_oa_glucose"), width=3),
                dbc.Col(dcc.Graph(id="at_oa_succinate+glucose"), width=3),
            ]
        ),
    ]
)


# --- Callback uses dictionary only ---
@app.callback(
    Output("at_succinate", "figure"),
    Output("oa_succinate", "figure"),
    Output("at_oa_succinate", "figure"),
    Output("at_glucose", "figure"),
    Output("oa_glucose", "figure"),
    Output("at_oa_glucose", "figure"),
    Output("at_succinate+glucose", "figure"),
    Output("oa_succinate+glucose", "figure"),
    Output("at_oa_succinate+glucose", "figure"),
    Input("μ_at_succ", "value"),
    Input("μ_at_gluc", "value"),
    Input("μ_oa_succ", "value"),
    Input("μ_oa_gluc", "value"),
    Input("at_K_succ", "value"),
    Input("at_K_gluc", "value"),
    Input("oa_K_succ", "value"),
    Input("oa_K_gluc", "value"),
    Input("at_a", "value"),
    Input("oa_a", "value"),
    Input("D", "value"),
)
def update_figures(
    μ_at_succ,
    μ_at_gluc,
    μ_oa_succ,
    μ_oa_gluc,
    at_K_succ,
    at_K_gluc,
    oa_K_succ,
    oa_K_gluc,
    at_a,
    oa_a,
    D,
):
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)
    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])
    at.v_succ = μ_at_succ
    at.v_gluc = μ_at_gluc
    oa.v_succ = μ_oa_succ
    oa.v_gluc = μ_oa_gluc
    at.K_succ = at_K_succ
    at.K_gluc = at_K_gluc
    oa.K_succ = oa_K_succ
    oa.K_gluc = oa_K_gluc
    at.a = at_a
    oa.a = oa_a
    D = D
    xs = np.linspace(0, 200, 500)

    # At mono-culture
    oa.N0 = 0.0
    at.N0 = 0.1
    at.a = 1
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_succ = model.plot_at_oa()

    at.a = 0
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_gluc = model.plot_at_oa()

    at.a = at_a
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_succ_gluc = model.plot_at_oa()

    # Oa mono-culture
    at.N0 = 0.0
    oa.N0 = 0.1
    oa.a = 1
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_oa_succ = model.plot_at_oa()

    oa.a = 0
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_oa_gluc = model.plot_at_oa()

    oa.a = oa_a
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_oa_succ_gluc = model.plot_at_oa()

    # At + Oa co-culture
    at.N0, oa.N0 = 0.05, 0.05
    at.a = 1
    oa.a = 1
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_oa_succ = model.plot_at_oa()

    at.a = 0
    oa.a = 0
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_oa_gluc = model.plot_at_oa()

    at.a = at_a
    oa.a = oa_a
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_oa_succ_gluc = model.plot_at_oa()

    return (
        fig_at_succ,
        fig_oa_succ,
        fig_at_oa_succ,
        fig_at_gluc,
        fig_oa_gluc,
        fig_at_oa_gluc,
        fig_at_succ_gluc,
        fig_oa_succ_gluc,
        fig_at_oa_succ_gluc,
    )


if __name__ == "__main__":
    app.run(debug=True)
