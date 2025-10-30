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


p_f = path.join("parameters", f"parameters_{conc}_mM_C.csv")
params = pd.read_csv(p_f, index_col=0)
at = Species("At", params.loc["At"])
oa = Species("Oa", params.loc["Oa"])


# --- Initialize parameters once ---
p_init = parse_params("parameters/parameters_dash_init.csv")


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
                dbc.Col(html.Span("v"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="p_v1",
                        type="number",
                        value=p_init["v1"],
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("v"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="p_v2",
                        type="number",
                        value=p_init["v2"],
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Span("K Resource A"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="p_K1_1",
                        type="number",
                        value=p_init["K1_1"],
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("K Resource A"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="p_K2_1",
                        type="number",
                        value=p_init["K2_1"],
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.Span("K Resource B"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="p_K1_2",
                        type="number",
                        value=p_init["K1_2"],
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("K Resource B"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="p_K2_2",
                        type="number",
                        value=p_init["K2_2"],
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
                        id="p_a1",
                        type="number",
                        value=p_init["a1"],
                        style={"width": "80px"},
                    ),
                    width=2,
                ),
                dbc.Col(html.Span("Resource preference"), width=2),
                dbc.Col(
                    dcc.Input(
                        id="p_a2",
                        type="number",
                        value=p_init["a2"],
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
    Output("at_succinate+glucose_chain", "figure"),
    Output("oa_succinate+glucose", "figure"),
    Output("oa_succinate+glucose_chain", "figure"),
    Output("at_oa_succinate+glucose", "figure"),
    Output("at_oa_succinate+glucose_chain", "figure"),
    Input("p_v1", "value"),
    Input("p_v2", "value"),
    Input("p_K1_1", "value"),
    Input("p_K2_1", "value"),
    Input("p_K1_2", "value"),
    Input("p_K2_2", "value"),
    Input("p_a1", "value"),
    Input("p_a2", "value"),
)
def update_figures(v1, v2, K1_1, K2_1, K1_2, K2_2, a1, a2):
    # start from the initial parameter dictionary
    base = p_init.copy()
    base.update(
        {
            "v1": v1,
            "v2": v2,
            "K1_1": K1_1,
            "K2_1": K2_1,
            "K1_2": K1_2,
            "K2_2": K2_2,
            "a1": a1,
            "a2": a2,
        }
    )

    # --- Condition 1: Succinate only ---


if __name__ == "__main__":
    app.run(debug=True)
