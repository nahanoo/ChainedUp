from dash import Dash, dcc, html, Output, Input
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from style import *
from plotly.subplots import make_subplots


def parse_params(f):
    df = dict(pd.read_csv(f))
    return pd.Series(df["value"].values, index=df["parameter"]).to_dict()


# --- Initialize parameters once ---
p_init = parse_params("parameters/parameters_dash_init.csv")


# --- Simple Monod-like model ---
def model(y, t, p):
    N1, N2, C1, C2, U1, U2, S1, S2 = y
    J1 = p["v1"] * (
        p["a1"] * C1 / (C1 + p["K1_1"]) + (1 - p["a1"]) * C2 / (C2 + p["K1_2"])
    )
    J2 = p["v2"] * (
        p["a2"] * C1 / (C1 + p["K2_1"]) + (1 - p["a2"]) * C2 / (C2 + p["K2_2"])
    )

    dN1 = J1 * N1 - p["D"] * N1
    dN2 = J2 * N2 - p["D"] * N2

    dC1 = (
        -p["v1"] * p["a1"] * C1 / (C1 + p["K1_1"]) * N1 / p["q1_1"]
        - p["v2"] * p["a2"] * C1 / (C1 + p["K2_1"]) * N2 / p["q2_1"]
        - p["D"] * C1
        + p["D"] * p["M1"]
    )
    dC2 = (
        -p["v1"] * (1 - p["a1"]) * C2 / (C2 + p["K1_2"]) * N1 / p["q1_2"]
        - p["v2"] * (1 - p["a2"]) * C2 / (C2 + p["K2_2"]) * N2 / p["q2_2"]
        - p["D"] * C2
        + p["D"] * p["M2"]
    )

    J1 = p["v1"] * (
        p["a1"] * S1 / (S1 + p["K1_1"]) + (1 - p["a1"]) * S2 / (S2 + p["K1_2"])
    )
    J2 = p["v2"] * (
        p["a2"] * S1 / (S1 + p["K2_1"]) + (1 - p["a2"]) * S2 / (S2 + p["K2_2"])
    )
    dU1 = J1 * U1 - p["D"] * U1 + p["D"] * N1
    dU2 = J2 * U2 - p["D"] * U2 + p["D"] * N2
    dS1 = (
        -p["v1"] * p["a1"] * S1 / (S1 + p["K1_1"]) * U1 / p["q1_1"]
        - p["v2"] * p["a2"] * S1 / (S1 + p["K2_1"]) * U2 / p["q2_1"]
        - p["D"] * S1
        + p["D"] * C1
    )
    dS2 = (
        -p["v1"] * (1 - p["a1"]) * S2 / (S2 + p["K1_2"]) * U1 / p["q1_2"]
        - p["v2"] * (1 - p["a2"]) * S2 / (S2 + p["K2_2"]) * U2 / p["q2_2"]
        - p["D"] * S2
        + p["D"] * C2
    )
    return [dN1, dN2, dC1, dC2, dU1, dU2, dS1, dS2]


def plot_species(y0, title, p):
    """Generate a growth curve using parameter dictionary p."""
    t = np.linspace(0, 100, 1000)
    N1, N2, C1, C2, U1, U2, S1, S2 = odeint(model, y0, t, args=(p,)).T

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Scatter(
            x=t,
            y=N1,
            mode="lines",
            name="Species At",
            legendgroup="Species At",
            line=dict(color=colors["at"], width=4),
            visible="legendonly" if np.max(N1) <= 1e-6 else True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t,
            y=N2,
            mode="lines",
            name="Species Oa",
            legendgroup="Species Oa",
            line=dict(color=colors["oa"], width=4),
            visible="legendonly" if np.max(N2) <= 1e-6 else True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=t,
            y=C1,
            mode="lines",
            name="Succinate",
            legendgroup="Succinate",
            line=dict(color="blue", width=4),
            visible="legendonly" if np.max(C1) <= 1e-6 else True,
        ),
        secondary_y=True,
    )
    fig.add_trace(
        go.Scatter(
            x=t,
            y=C2,
            mode="lines",
            name="Glucose",
            legendgroup="Glucose",
            line=dict(color="orange", width=4),
            visible="legendonly" if np.max(C2) <= 1e-6 else True,
        ),
        secondary_y=True,
    )
    fig.update_layout(
        title=title,
        height=300,
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD600", ticks="inside", range=[0.0, 1]),
        yaxis2=dict(title="Resource concentration", ticks="inside", range=[0.0, 1]),
    )
    fig = style_plot(fig)

    fig_chain = make_subplots(specs=[[{"secondary_y": True}]])
    fig_chain.add_trace(
        go.Scatter(
            x=t,
            y=N1,
            mode="lines",
            name="Species At",
            legendgroup="Species At",
            line=dict(color=colors["at"], width=4),
            visible="legendonly" if np.max(N1) <= 1e-6 else True,
        )
    )
    fig_chain.add_trace(
        go.Scatter(
            x=t,
            y=N2,
            mode="lines",
            name="Species Oa",
            legendgroup="Species Oa",
            line=dict(color=colors["oa"], width=4),
            visible="legendonly" if np.max(N2) <= 1e-6 else True,
        )
    )
    fig_chain.add_trace(
        go.Scatter(
            x=t,
            y=C1,
            mode="lines",
            name="Succinate",
            legendgroup="Succinate",
            line=dict(color="blue", width=4),
            visible="legendonly" if np.max(C1) <= 1e-6 else True,
        ),
        secondary_y=True,
    )
    fig_chain.add_trace(
        go.Scatter(
            x=t,
            y=C2,
            mode="lines",
            name="Glucose",
            legendgroup="Glucose",
            line=dict(color="orange", width=4),
            visible="legendonly" if np.max(C2) <= 1e-6 else True,
        ),
        secondary_y=True,
    )
    fig_chain.update_layout(
        title=title,
        height=300,
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD600", ticks="inside", range=[0.0, 1]),
        yaxis2=dict(title="Resource concentration", ticks="inside", range=[0.0, 1]),
    )

    fig_chain = style_plot(fig_chain)
    return fig, fig_chain


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
                dbc.Col(dcc.Graph(id="at_succinate+glucose_chain"), width=3),
            ]
        ),
        dbc.Row(html.H4("Oa mono-culture")),
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id="oa_succinate"), width=3),
                dbc.Col(dcc.Graph(id="oa_glucose"), width=3),
                dbc.Col(dcc.Graph(id="oa_succinate+glucose"), width=3),
                dbc.Col(dcc.Graph(id="oa_succinate+glucose_chain"), width=3),
            ]
        ),
        dbc.Row(html.H4("At + Oa co-culture")),
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id="at_oa_succinate"), width=3),
                dbc.Col(dcc.Graph(id="at_oa_glucose"), width=3),
                dbc.Col(dcc.Graph(id="at_oa_succinate+glucose"), width=3),
                dbc.Col(dcc.Graph(id="at_oa_succinate+glucose_chain"), width=3),
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
    p_succ = base.copy()
    p_succ["M1"], p_succ["M2"] = 1, 0
    p_succ["a1"] = 1
    p_succ["a2"] = 1
    fig1 = plot_species([0.1, 0, 1, 0, 0, 0, 0, 0], "At Succinate", p_succ)[0]
    fig1.update_layout(showlegend=False)

    fig2 = plot_species([0, 0.1, 1, 0, 0, 0, 0, 0], "Oa Succinate", p_succ)[0]
    fig2.update_layout(showlegend=False)

    fig3 = plot_species([0.1, 0.1, 1, 0, 0, 0, 0, 0], "At + Oa Succinate", p_succ)[0]
    fig3.update_layout(showlegend=False)

    # --- Condition 2: Glucose only ---
    p_glu = base.copy()
    p_glu["M1"], p_glu["M2"] = 0, 1
    p_glu["a1"] = 0
    p_glu["a2"] = 0
    fig4 = plot_species([0.1, 0, 0, 1, 0, 0, 0, 0], "At Glucose", p_glu)[0]
    fig4.update_layout(showlegend=False)
    fig5 = plot_species([0, 0.1, 0, 1, 0, 0, 0, 0], "Oa Glucose", p_glu)[0]
    fig5.update_layout(showlegend=False)
    fig6 = plot_species([0.1, 0.1, 0, 1, 0, 0, 0, 0], "At + Oa Glucose", p_glu)[0]
    fig6.update_layout(showlegend=False)

    # --- Condition 3: Both substrates ---
    p_mix = base.copy()
    p_mix["M1"], p_mix["M2"] = 1, 1
    fig7 = plot_species([0.1, 0, 1, 1, 0, 0, 0, 0], "At Succinate + Glucose", p_mix)[0]
    fig7.update_layout(showlegend=False)
    fig8 = plot_species([0.1, 0, 1, 1, 0, 0, 0, 0], "At Succinate + Glucose", p_mix)[1]
    fig9 = plot_species([0, 0.1, 1, 1, 0, 0, 0, 0], "Oa Succinate + Glucose", p_mix)[0]
    fig9.update_layout(showlegend=False)
    fig10 = plot_species([0, 0.1, 1, 1, 0, 0, 0, 0], "Oa Succinate + Glucose", p_mix)[1]
    fig11 = plot_species(
        [0.1, 0.1, 1, 1, 0, 0, 0, 0], "At + Oa Succinate + Glucose", p_mix
    )[0]
    fig11.update_layout(showlegend=False)
    fig12 = plot_species(
        [0.1, 0.1, 1, 1, 0, 0, 0, 0], "At + Oa Succinate + Glucose", p_mix
    )[1]

    return fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12


if __name__ == "__main__":
    app.run(debug=True)
