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

# --------------------
# Initial setup
# --------------------
conc = 15
p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)
at = Species("At", params.loc["At"])
oa = Species("Oa", params.loc["Oa"])
global D
D = 0.1
M = Model(at, oa, None, np.linspace(0, 200, 500), conc, D)

# --------------------
# Dash app
# --------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])


# Helper for resource figures: robust even if Model has no plot_resources()
def get_resource_fig(model, title: str) -> go.Figure:
    if hasattr(model, "plot_resources"):
        fig = model.plot_resources()
    else:
        # Fallback: avoid crashing if plot_resources() is not implemented
        fig = go.Figure()
    fig.update_layout(title=title)
    return fig


# --------------------
# Layout
# --------------------
app.layout = dbc.Container(
    [
        dbc.Row(html.H3("Parameters")),
        dbc.Row(
            [
                dbc.Col([html.H4("At")], width=4),
                dbc.Col(html.H4("Oa"), width=4),
            ]
        ),
        # μ controls
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
        # K controls
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
        # a controls
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
        # D control
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
        # ---- TABS ----
        dbc.Tabs(
            id="main-tabs",
            active_tab="tab-biomass",
            children=[
                # Tab 1: Biomass
                dbc.Tab(
                    label="Biomass",
                    tab_id="tab-biomass",
                    children=[
                        dbc.Row(html.H4("At mono-culture")),
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="at_succinate"), width=3),
                                dbc.Col(dcc.Graph(id="at_glucose"), width=3),
                                dbc.Col(
                                    dcc.Graph(id="at_succinate+glucose"),
                                    width=3,
                                ),
                            ]
                        ),
                        dbc.Row(html.H4("Oa mono-culture")),
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="oa_succinate"), width=3),
                                dbc.Col(dcc.Graph(id="oa_glucose"), width=3),
                                dbc.Col(
                                    dcc.Graph(id="oa_succinate+glucose"),
                                    width=3,
                                ),
                            ]
                        ),
                        dbc.Row(html.H4("At + Oa co-culture")),
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="at_oa_succinate"), width=3),
                                dbc.Col(dcc.Graph(id="at_oa_glucose"), width=3),
                                dbc.Col(
                                    dcc.Graph(id="at_oa_succinate+glucose"),
                                    width=3,
                                ),
                            ]
                        ),
                    ],
                ),
                # Tab 2: Resources
                dbc.Tab(
                    label="Resources",
                    tab_id="tab-resources",
                    children=[
                        dbc.Row(html.H4("At mono-culture (resources)")),
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="at_succinate_res"), width=3),
                                dbc.Col(dcc.Graph(id="at_glucose_res"), width=3),
                                dbc.Col(
                                    dcc.Graph(id="at_succinate+glucose_res"),
                                    width=3,
                                ),
                            ]
                        ),
                        dbc.Row(html.H4("Oa mono-culture (resources)")),
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="oa_succinate_res"), width=3),
                                dbc.Col(dcc.Graph(id="oa_glucose_res"), width=3),
                                dbc.Col(
                                    dcc.Graph(id="oa_succinate+glucose_res"),
                                    width=3,
                                ),
                            ]
                        ),
                        dbc.Row(html.H4("At + Oa co-culture (resources)")),
                        dbc.Row(
                            [
                                dbc.Col(dcc.Graph(id="at_oa_succinate_res"), width=3),
                                dbc.Col(dcc.Graph(id="at_oa_glucose_res"), width=3),
                                dbc.Col(
                                    dcc.Graph(id="at_oa_succinate+glucose_res"),
                                    width=3,
                                ),
                            ]
                        ),
                    ],
                ),
            ],
        ),
    ]
)


# --------------------
# Callback
# --------------------
@app.callback(
    # Biomass plots
    Output("at_succinate", "figure"),
    Output("oa_succinate", "figure"),
    Output("at_oa_succinate", "figure"),
    Output("at_glucose", "figure"),
    Output("oa_glucose", "figure"),
    Output("at_oa_glucose", "figure"),
    Output("at_succinate+glucose", "figure"),
    Output("oa_succinate+glucose", "figure"),
    Output("at_oa_succinate+glucose", "figure"),
    # Resource plots
    Output("at_succinate_res", "figure"),
    Output("oa_succinate_res", "figure"),
    Output("at_oa_succinate_res", "figure"),
    Output("at_glucose_res", "figure"),
    Output("oa_glucose_res", "figure"),
    Output("at_oa_glucose_res", "figure"),
    Output("at_succinate+glucose_res", "figure"),
    Output("oa_succinate+glucose_res", "figure"),
    Output("at_oa_succinate+glucose_res", "figure"),
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
    D_in,
):
    conc = 15
    p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
    params = pd.read_csv(p_f, index_col=0)

    at = Species("At", params.loc["At"])
    oa = Species("Oa", params.loc["Oa"])

    # Update parameters from inputs
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

    xs = np.linspace(0, 200, 500)
    D = 0.0 if D_in is None else D_in

    # --------------------
    # At mono-culture
    # --------------------
    oa.N0 = 0.0
    at.N0 = 0.1

    # At on succinate
    at.a = 1
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_succ = model.plot_at_oa()
    fig_at_succ.update_layout(
        title=f"At on {model.C_to_mM_succinate[model.C_mono]} mM Succinate"
    )
    fig_at_succ_res = get_resource_fig(
        model,
        title=(
            f"Resources: At on {model.C_to_mM_succinate[model.C_mono]} mM Succinate"
        ),
    )

    # At on glucose
    at.a = 0
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_gluc = model.plot_at_oa()
    fig_at_gluc.update_layout(title="At on Glucose")
    fig_at_gluc_res = get_resource_fig(model, title="Resources: At on Glucose")

    # At on succinate + glucose
    at.a = at_a
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_succ_gluc = model.plot_at_oa()
    fig_at_succ_gluc.update_layout(
        title=(
            f"At on {model.C_to_mM_succinate[model.C_mono]} mM Succinate + "
            f"{model.C_to_mM_glucose[model.C_mono]} mM Glucose"
        )
    )
    fig_at_succ_gluc_res = get_resource_fig(
        model,
        title=(
            f"Resources: At on {model.C_to_mM_succinate[model.C_mono]} mM Succinate + "
            f"{model.C_to_mM_glucose[model.C_mono]} mM Glucose"
        ),
    )

    # --------------------
    # Oa mono-culture
    # --------------------
    at.N0 = 0.0
    oa.N0 = 0.1

    # Oa on succinate
    oa.a = 1
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_oa_succ = model.plot_at_oa()
    fig_oa_succ.update_layout(
        title=f"Oa on {model.C_to_mM_succinate[model.C_mono]} mM Succinate"
    )
    fig_oa_succ_res = get_resource_fig(
        model,
        title=(
            f"Resources: Oa on {model.C_to_mM_succinate[model.C_mono]} mM Succinate"
        ),
    )

    # Oa on glucose
    oa.a = 0
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_oa_gluc = model.plot_at_oa()
    fig_oa_gluc.update_layout(
        title=f"Oa on {model.C_to_mM_glucose[model.C_mono]} mM Glucose"
    )
    fig_oa_gluc_res = get_resource_fig(
        model,
        title=(f"Resources: Oa on {model.C_to_mM_glucose[model.C_mono]} mM Glucose"),
    )

    # Oa on succinate + glucose
    oa.a = oa_a
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_oa_succ_gluc = model.plot_at_oa()
    fig_oa_succ_gluc.update_layout(
        title=(
            f"Oa on {model.C_to_mM_succinate[model.C_mono]} mM Succinate + "
            f"{model.C_to_mM_glucose[model.C_mono]} mM Glucose"
        )
    )
    fig_oa_succ_gluc_res = get_resource_fig(
        model,
        title=(
            f"Resources: Oa on {model.C_to_mM_succinate[model.C_mono]} mM Succinate + "
            f"{model.C_to_mM_glucose[model.C_mono]} mM Glucose"
        ),
    )

    # --------------------
    # At + Oa co-culture
    # --------------------
    at.N0, oa.N0 = 0.05, 0.05

    # Co-culture on succinate
    at.a = 1
    oa.a = 1
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_oa_succ = model.plot_at_oa()
    fig_at_oa_succ.update_layout(
        title=f"At + Oa on {model.C_to_mM_succinate[model.C_mono]} mM Succinate"
    )
    fig_at_oa_succ_res = get_resource_fig(
        model,
        title=(
            f"Resources: At + Oa on "
            f"{model.C_to_mM_succinate[model.C_mono]} mM Succinate"
        ),
    )

    # Co-culture on glucose
    at.a = 0
    oa.a = 0
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_oa_gluc = model.plot_at_oa()
    fig_at_oa_gluc.update_layout(
        title=f"At + Oa on {model.C_to_mM_glucose[model.C_mono]} mM Glucose"
    )
    fig_at_oa_gluc_res = get_resource_fig(
        model,
        title=(
            f"Resources: At + Oa on "
            f"{model.C_to_mM_glucose[model.C_mono]} mM Glucose"
        ),
    )

    # Co-culture on succinate + glucose
    at.a = at_a
    oa.a = oa_a
    model = Model(at, oa, None, xs, conc, D)
    model.integrate_model()
    fig_at_oa_succ_gluc = model.plot_at_oa()
    fig_at_oa_succ_gluc.update_layout(
        title=(
            f"At + Oa on {model.C_to_mM_succinate[model.C_mono]} mM Succinate + "
            f"{model.C_to_mM_glucose[model.C_mono]} mM Glucose"
        )
    )
    fig_at_oa_succ_gluc_res = get_resource_fig(
        model,
        title=(
            f"Resources: At + Oa on "
            f"{model.C_to_mM_succinate[model.C_mono]} mM Succinate + "
            f"{model.C_to_mM_glucose[model.C_mono]} mM Glucose"
        ),
    )

    return (
        # biomass
        fig_at_succ,
        fig_oa_succ,
        fig_at_oa_succ,
        fig_at_gluc,
        fig_oa_gluc,
        fig_at_oa_gluc,
        fig_at_succ_gluc,
        fig_oa_succ_gluc,
        fig_at_oa_succ_gluc,
        # resources
        fig_at_succ_res,
        fig_oa_succ_res,
        fig_at_oa_succ_res,
        fig_at_gluc_res,
        fig_oa_gluc_res,
        fig_at_oa_glucose_res,
        fig_at_succ_gluc_res,
        fig_oa_succ_gluc_res,
        fig_at_oa_succ_gluc_res,
    )


if __name__ == "__main__":
    app.run(debug=True)
