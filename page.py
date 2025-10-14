from dash import Dash, dcc, html, Output, Input
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from style import *


def parse_params(cfus=False):
    if cfus:
        df = dict(pd.read_csv("parameters_cfus.csv"))
        params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
        p = params
        return p
    else:
        df = dict(pd.read_csv("parameters.csv"))
        params = pd.Series(df["value"].values, index=df["parameter"]).to_dict()
        p = params
        return p


# --- Simple Monod-like model ---
def model(y, t, p, a, b):
    N1, N2, C1, C2 = y
    J1 = p["v1_1"] * (a * C1 / (C1 + p["K1_1"]) + (1 - a) * C2 / (C2 + p["K1_2"]))
    J2 = p["v2_1"] * (b * C1 / (C1 + p["K2_1"]) + (1 - b) * C2 / (C2 + p["K2_2"]))
    """if J1 < p["D"]:
        J1 = p["D"]
    if J2 < p["D"]:
        J2 = p["D"]"""

    dN1 = J1 * N1 - p["D"] * N1
    dN2 = J2 * N2 - p["D"] * N2
    dC1 = (
        -p["v1_1"] * a * C1 / (C1 + p["K1_1"]) * N1 / p["q1_1"]
        - p["v2_1"] * b * C1 / (C1 + p["K2_1"]) * N2 / p["q2_1"]
        - p["D"] * C1
        + p["D"] * p["M1"]
    )
    dC2 = (
        -p["v1_1"] * (1 - a) * C2 / (C2 + p["K1_2"]) * N1 / p["q1_2"]
        - p["v2_1"] * (1 - b) * C2 / (C2 + p["K2_2"]) * N2 / p["q2_2"]
        - p["D"] * C2
        + p["D"] * p["M2"]
    )

    return [dN1, dN2, dC1, dC2]


xs = np.linspace(0, 200, 1000)


def mono_culture_single_resource(N1, N2, M1, M2, a, b, title):
    p = parse_params(cfus=True)
    p["N1"] = N1
    p["N2"] = N2
    p["M1"] = M1
    p["M2"] = M2
    N1, N2, C1, C2 = odeint(
        model,
        [p["N1"], p["N2"], p["M1"], p["M2"]],
        xs,
        args=(p, a, b),
    ).T

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=N1,
            mode="lines",
            line=dict(width=3, color=colors["at"]),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=N2,
            mode="lines",
            line=dict(width=3, color=colors["oa"]),
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        title=title,
        yaxis=dict(
            title="Population density [CFU/mL]",
            ticks="inside",
            range=[4, 9],
            dtick=0.2,
            type="log",
            tickformat=".0e",
        ),
    )
    fig = style_plot(fig, line_thickness=4, font_size=20)

    return fig


# --- Initialize app ---
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# --- Layout ---
app.layout = dbc.Container(
    [
        html.H3("Resource preferences", className="mt-3 mb-4"),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label(
                            "Species 1 (1 = only resource A, 0 = only resource B)",
                            className="mt-3",
                            style={"font-size": "18px", "color": "black"},
                        ),
                        dcc.Slider(
                            id="slider-a",
                            min=0,
                            max=1,
                            step=0.1,
                            value=1.0,
                            marks={
                                v: {
                                    "label": f"{v:.1f}",
                                    "style": {"color": "black", "font-size": "16px"},
                                }
                                for v in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
                            },
                            tooltip={
                                "always_visible": True,
                                "placement": "bottom",
                                "style": {
                                    "color": "black",
                                    "fontSize": "16px",
                                },
                            },
                        ),
                    ],
                    width=2,
                ),
                dbc.Col(
                    [
                        html.Label(
                            "Species 2 (1 = only resource A, 0 = only resource B)",
                            className="mt-3",
                            style={"font-size": "18px", "color": "black"},
                        ),
                        dcc.Slider(
                            id="slider-b",
                            min=0,
                            max=1,
                            step=0.1,
                            value=1.0,
                            marks={
                                v: {
                                    "label": f"{v:.1f}",
                                    "style": {"color": "black", "font-size": "16px"},
                                }
                                for v in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
                            },
                            tooltip={
                                "always_visible": True,
                                "placement": "bottom",
                                "style": {
                                    "color": "black",
                                    "fontSize": "16px",
                                },
                            },
                        ),
                    ],
                    width=2,
                ),
            ]
        ),
        html.H3("Species 1", className="mt-3 mb-4"),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Graph(
                        id="resource-a",
                        figure=mono_culture_single_resource(
                            1e7, 0, 1, 0, 1, 1, "Resource A"
                        ),
                        style={"height": "30vh"},
                    ),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(
                        id="resource-b",
                        figure=mono_culture_single_resource(
                            1e7, 0, 0, 1, 0, 1, "Resource B"
                        ),
                        style={"height": "30vh"},
                    ),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(id="resource-a-b", style={"height": "30vh"}),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(id="resource-conc-a-b", style={"height": "30vh"}),
                    width=3,
                ),
            ],
        ),
        html.H3("Species 2", className="mt-3 mb-4"),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Graph(
                        id="resource-a-species-2",
                        figure=mono_culture_single_resource(
                            0, 1e7, 1, 0, 1, 1, "Resource A"
                        ),
                        style={"height": "30vh"},
                    ),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(
                        id="resource-b-species-2",
                        figure=mono_culture_single_resource(
                            0, 1e7, 0, 1, 1, 0, "Resource B"
                        ),
                        style={"height": "30vh"},
                    ),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(id="resource-a-b-species-2", style={"height": "30vh"}),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(
                        id="resource-conc-a-b-species-2", style={"height": "30vh"}
                    ),
                    width=3,
                ),
            ],
        ),
        html.H3("Species 1 + 2", className="mt-3 mb-4"),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Graph(
                        id="resource-a-species-1-2",
                        figure=mono_culture_single_resource(
                            1e7, 1e7, 1, 0, 1, 1, "Resource A"
                        ),
                        style={"height": "30vh"},
                    ),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(
                        id="resource-b-species-1-2",
                        figure=mono_culture_single_resource(
                            1e7, 1e7, 0, 1, 0, 0, "Resource B"
                        ),
                        style={"height": "30vh"},
                    ),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(id="resource-a-b-species-1-2", style={"height": "30vh"}),
                    width=3,
                ),
                dbc.Col(
                    dcc.Graph(
                        id="resource-conc-a-b-species-1-2", style={"height": "30vh"}
                    ),
                    width=3,
                ),
            ],
        ),
    ],
    fluid=True,
)


@app.callback(
    Output("resource-a-b", "figure"),
    Output("resource-conc-a-b", "figure"),
    Input("slider-a", "value"),
    Input("slider-b", "value"),
)
def update_plot(a, b):
    p = parse_params(cfus=True)
    p["N2"] = 0
    N1, N2, C1, C2 = odeint(
        model,
        [p["N1"], p["N2"], p["M1"], p["M2"]],
        xs,
        args=(p, a, b),
    ).T

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=N1,
            mode="lines",
            line=dict(width=3, color=colors["at"]),
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Population density [CFU/mL]",
            ticks="inside",
            range=[4, 9],
            dtick=0.2,
            type="log",
            tickformat=".0e",
        ),
        title="Resource A + B",
    )
    fig = style_plot(fig, line_thickness=3, font_size=16)

    fig_resource = go.Figure()
    fig_resource.add_trace(
        go.Scatter(
            x=xs, y=C1, mode="lines", line=dict(width=3, color="blue"), name="C1"
        )
    )
    fig_resource.add_trace(
        go.Scatter(
            x=xs, y=C2, mode="lines", line=dict(width=3, color="orange"), name="C2"
        )
    )

    fig_resource.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Resource concentration",
            ticks="inside",
            # range=[7, 9],
            # dtick=1,
            # type="log",
            # tickformat=".0e",
        ),
    )
    fig_resource = style_plot(fig_resource, line_thickness=3, font_size=16)

    return fig, fig_resource


@app.callback(
    Output("resource-a-b-species-2", "figure"),
    Output("resource-conc-a-b-species-2", "figure"),
    Input("slider-a", "value"),
    Input("slider-b", "value"),
)
def update_plot(a, b):
    p = parse_params(cfus=True)
    p["N1"] = 0
    N1, N2, C1, C2 = odeint(
        model,
        [p["N1"], p["N2"], p["M1"], p["M2"]],
        xs,
        args=(p, a, b),
    ).T

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=N2,
            mode="lines",
            line=dict(width=3, color=colors["oa"]),
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Population density [CFU/mL]",
            ticks="inside",
            range=[4, 9],
            dtick=0.2,
            type="log",
            tickformat=".0e",
        ),
        title="Resource A + B",
    )
    fig = style_plot(fig, line_thickness=3, font_size=16)

    fig_resource = go.Figure()
    fig_resource.add_trace(
        go.Scatter(
            x=xs, y=C1, mode="lines", line=dict(width=3, color="blue"), name="C1"
        )
    )
    fig_resource.add_trace(
        go.Scatter(
            x=xs, y=C2, mode="lines", line=dict(width=3, color="orange"), name="C2"
        )
    )

    fig_resource.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Resource concentration",
            ticks="inside",
            # range=[7, 9],
            # dtick=1,
            # type="log",
            # tickformat=".0e",
        ),
    )
    fig_resource = style_plot(fig_resource, line_thickness=3, font_size=16)

    return fig, fig_resource


@app.callback(
    Output("resource-a-b-species-1-2", "figure"),
    Output("resource-conc-a-b-species-1-2", "figure"),
    Input("slider-a", "value"),
    Input("slider-b", "value"),
)
def update_plot(a, b):
    p = parse_params(cfus=True)
    N1, N2, C1, C2 = odeint(
        model,
        [p["N1"], p["N2"], p["M1"], p["M2"]],
        xs,
        args=(p, a, b),
    ).T

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=N1,
            mode="lines",
            line=dict(width=3, color=colors["at"]),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=N2,
            mode="lines",
            line=dict(width=3, color=colors["oa"]),
        )
    )

    fig.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Population density [CFU/mL]",
            ticks="inside",
            range=[4, 9],
            dtick=0.2,
            type="log",
            tickformat=".0e",
        ),
        title="Resource A + B",
    )
    fig = style_plot(fig, line_thickness=3, font_size=16)

    fig_resource = go.Figure()
    fig_resource.add_trace(
        go.Scatter(
            x=xs, y=C1, mode="lines", line=dict(width=3, color="blue"), name="C1"
        )
    )
    fig_resource.add_trace(
        go.Scatter(
            x=xs, y=C2, mode="lines", line=dict(width=3, color="orange"), name="C2"
        )
    )

    fig_resource.update_layout(
        xaxis=dict(
            title="Time [h]",
            ticks="inside",
        ),
        yaxis=dict(
            title="Resource concentration",
            ticks="inside",
            # range=[7, 9],
            # dtick=1,
            # type="log",
            # tickformat=".0e",
        ),
    )
    fig_resource = style_plot(fig_resource, line_thickness=3, font_size=16)

    return fig, fig_resource


# --- Run ---
if __name__ == "__main__":
    app.run(debug=True)
