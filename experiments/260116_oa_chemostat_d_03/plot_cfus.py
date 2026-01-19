from style import colors, style_plot
from chibio_parser import cfu_parser
import plotly.graph_objects as go


def plot_all_conditions():
    df, _ = cfu_parser("251208_at_chemostat_d_03")
    df["average"] = df["average"].replace(0, None)
    reactor_map = {
        "M0": "Succinate",
        "M1": "Glucose",
        "M2": "Succinate+Glucose",
        "M3": "Succinate+Glucose Outflow",
    }

    fig = go.Figure()
    reactors = sorted(df["reactor"].unique())

    for r in reactors:
        medium = reactor_map[r]
        dfr = df[df["reactor"] == r].sort_values("sample_time")

        fig.add_trace(
            go.Scatter(
                x=dfr["sample_time"],
                y=dfr["average"],
                error_y=dict(type="data", array=dfr["stdev"], visible=True),
                mode="lines+markers",
                name=medium,
                marker=dict(color=colors[medium], size=6),
                line=dict(color=colors[medium]),
            )
        )

    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        yaxis=dict(ticks="inside", title="CFU/mL", type="log", range=[6, 10]),
        title="CFU/mL in chemostats at dilution rate 0.3 h<sup>-1</sup>",
        legend=dict(title="Carbon Source", x=0.00, y=0.01),
    )

    fig = style_plot(fig, line_thickness=2, marker_size=6, font_size=12)
    fig.write_image("plots/cfu_chemostat.svg")


def plot_succinate_glucose():
    df, _ = cfu_parser("251208_at_chemostat_d_03")
    df = df[(df["reactor"] == "M0") | (df["reactor"] == "M1")]
    df["average"] = df["average"].replace(0, None)
    reactor_map = {
        "M0": "Succinate",
        "M1": "Glucose",
        "M2": "Succinate+Glucose",
        "M3": "Succinate+Glucose Outflow",
    }

    fig = go.Figure()
    reactors = sorted(df["reactor"].unique())

    for r in reactors:
        medium = reactor_map[r]
        dfr = df[df["reactor"] == r].sort_values("sample_time")

        fig.add_trace(
            go.Scatter(
                x=dfr["sample_time"],
                y=dfr["average"],
                error_y=dict(type="data", array=dfr["stdev"], visible=True),
                mode="lines+markers",
                name=medium,
                marker=dict(color=colors[medium], size=6),
                line=dict(color=colors[medium]),
            )
        )

    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        yaxis=dict(ticks="inside", title="CFU/mL", type="log", range=[6, 10]),
        title="CFU/mL in chemostats at dilution rate 0.3 h<sup>-1</sup>",
        legend=dict(title="Carbon Source", x=0.00, y=0.01),
        height=190,
        width=540,
        showlegend=False,
    )

    fig = style_plot(fig, line_thickness=2, marker_size=6, font_size=14)
    fig.write_image("plots/cfu_chemostat_succinate_glucose.svg")


def plot_chain():
    df, _ = cfu_parser("251208_at_chemostat_d_03")
    df = df[(df["reactor"] == "M2") | (df["reactor"] == "M3")]
    df["average"] = df["average"].replace(0, None)
    reactor_map = {
        "M0": "Succinate",
        "M1": "Glucose",
        "M2": "Succinate+Glucose",
        "M3": "Succinate+Glucose Outflow",
    }

    fig = go.Figure()
    reactors = sorted(df["reactor"].unique())
    names = ["Chemostat 1", "Downstream Chemostat"]
    for r in reactors:
        medium = reactor_map[r]
        dfr = df[df["reactor"] == r].sort_values("sample_time")

        fig.add_trace(
            go.Scatter(
                x=dfr["sample_time"],
                y=dfr["average"],
                error_y=dict(type="data", array=dfr["stdev"], visible=True),
                mode="lines+markers",
                name=names[reactors.index(r)],
                marker=dict(color=colors[medium], size=6),
                line=dict(color=colors[medium]),
            )
        )

    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        yaxis=dict(ticks="inside", title="CFU/mL", type="log", range=[6, 10]),
        title="Mono-culture Succinate+Glucose",
        legend=dict(x=0.00, y=0.01),
        height=190,
        width=250,
        showlegend=False,
    )

    fig = style_plot(fig, line_thickness=2, marker_size=6, font_size=14)
    fig.write_image("plots/cfu_chemostat_succinate_glucose_chain.svg")


def plot_succinate_glucose_single():
    df, _ = cfu_parser("251208_at_chemostat_d_03")
    df = df[(df["reactor"] == "M2")]
    df["average"] = df["average"].replace(0, None)
    reactor_map = {
        "M0": "Succinate",
        "M1": "Glucose",
        "M2": "Succinate+Glucose",
        "M3": "Succinate+Glucose Outflow",
    }

    fig = go.Figure()
    reactors = sorted(df["reactor"].unique())

    for r in reactors:
        medium = reactor_map[r]
        dfr = df[df["reactor"] == r].sort_values("sample_time")

        fig.add_trace(
            go.Scatter(
                x=dfr["sample_time"],
                y=dfr["average"],
                error_y=dict(type="data", array=dfr["stdev"], visible=True),
                mode="lines+markers",
                name=medium,
                marker=dict(color=colors[medium], size=6),
                line=dict(color=colors[medium]),
            )
        )

    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        yaxis=dict(ticks="inside", title="CFU/mL", type="log", range=[6, 10]),
        title="Mono-culture Succinate+Glucose",
        legend=dict(title="Carbon Source", x=0.00, y=0.01),
        height=190,
        width=250,
        showlegend=False,
    )

    fig = style_plot(fig, line_thickness=2, marker_size=6, font_size=14)
    fig.write_image("plots/cfu_chemostat_succinate+glucose.svg")


plot_all_conditions()
plot_succinate_glucose()
plot_chain()
plot_succinate_glucose_single()
