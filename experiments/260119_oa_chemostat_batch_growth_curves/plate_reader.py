from experiment import Species, Experiment
import plotly.graph_objects as go
from style import colors, style_plot


def plot_spent_chemostat_growth_curves():
    def filter_lgs(lgs):
        filtered_lgs = []
        css = [
            "Succinate spent medium",
            "Glucose spent medium",
            "Succinate+Glucose spent medium",
            "Succinate+Glucose outflow spent medium",
        ]
        for lg in lgs:
            if lg.carbon_source in css:
                filtered_lgs.append(lg)
        return filtered_lgs

    e = Experiment(d="../../data/260119_oa_chemostat_batch_growth/metaod/")
    e.build_conditions()

    adjust_od = {
        "Succinate spent medium": 0.42,
        "Glucose spent medium": 0.15,
        "Succinate+Glucose spent medium": 0.75,
        "Succinate+Glucose outflow spent medium": 0.82,
    }

    color_map = {
        "Succinate spent medium": colors["Succinate"],
        "Glucose spent medium": colors["Glucose"],
        "Succinate+Glucose spent medium": colors["Succinate+Glucose"],
        "Succinate+Glucose outflow spent medium": colors["Succinate+Glucose Outflow"],
    }
    lgs = filter_lgs(e.conditions)
    fig = go.Figure()
    for lg in lgs:
        for x, y in zip(lg.xs, lg.ys):
            y0 = adjust_od[lg.carbon_source] - y[0]
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y + y0,
                    mode="lines",
                    name=f"{lg.carbon_source}",
                    showlegend=False,
                    line=dict(color=color_map[lg.carbon_source]),
                )
            )
    for color in color_map.keys():
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines+markers",
                name=f"{color}",
                line=dict(color=color_map[color]),
                marker=dict(size=6),
            )
        )
    fig.update_layout(
        title="Growth Curves following chemostat experiment",
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD600", ticks="inside"),
        legend=dict(
            title="Carbon Source",
            x=0.99,
            y=0.01,
            xanchor="right",
            yanchor="bottom",
            bgcolor="rgba(0,0,0,0)",
        ),
    )
    fig = style_plot(fig, marker_size=6, line_thickness=2.5, font_size=12)
    fig.write_image("plots/chemostat_growth_curve.svg")


def plot_fresh_chemostat_growth_curves():
    def filter_lgs(lgs):
        filtered_lgs = []
        css = [
            "Succinate fresh medium",
            "Glucose fresh medium",
            "Succinate+Glucose fresh medium",
            "Succinate+Glucose outflow fresh medium",
        ]
        for lg in lgs:
            if lg.carbon_source in css:
                filtered_lgs.append(lg)
        return filtered_lgs

    e = Experiment(d="../../data/260119_oa_chemostat_batch_growth/metaod/")
    e.build_conditions()

    adjust_od = {
        "Succinate fresh medium": 0.42,
        "Glucose fresh medium": 0.15,
        "Succinate+Glucose fresh medium": 0.75,
        "Succinate+Glucose outflow fresh medium": 0.82,
    }

    color_map = {
        "Succinate fresh medium": colors["Succinate"],
        "Glucose fresh medium": colors["Glucose"],
        "Succinate+Glucose fresh medium": colors["Succinate+Glucose"],
        "Succinate+Glucose outflow fresh medium": colors["Succinate+Glucose Outflow"],
    }
    lgs = filter_lgs(e.conditions)
    fig = go.Figure()
    for lg in lgs:
        for x, y in zip(lg.xs, lg.ys):
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    name=f"{lg.carbon_source}",
                    showlegend=False,
                    line=dict(color=color_map[lg.carbon_source]),
                )
            )
    for color in color_map.keys():
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines+markers",
                name=f"{color}",
                line=dict(color=color_map[color]),
                marker=dict(size=6),
            )
        )
    fig.update_layout(
        title="Growth Curves following chemostat experiment",
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD600", ticks="inside"),
        legend=dict(
            title="Carbon Source",
            x=0.99,
            y=0.01,
            xanchor="right",
            yanchor="bottom",
            bgcolor="rgba(0,0,0,0)",
        ),
    )
    fig = style_plot(fig, marker_size=6, line_thickness=2.5, font_size=12)
    fig.write_image("plots/chemostat_fresh_mediumgrowth_curve.svg")


def plot_growth_curves_stocks():
    fig = go.Figure()
    e = Experiment(d="../../data/251018_succinate_glucose_plate_reader/metaod/")
    e.build_conditions()
    carbon_sources = ["Succinate", "Glucose", "Succinate+Glucose"]
    c_to_mM = {
        5: ["1.25 mM Succinate", "0.833 mM Glucose"],
        15: ["3.75 mM Succinate", "2.5 mM Glucose"],
        30: ["7.5 mM Succinate", "5 mM Glucose"],
    }
    conc_init = 15
    concentrations = [conc_init, conc_init, 2 * conc_init]
    for cs, conc in zip(carbon_sources, concentrations):
        for c in e.conditions:
            if (
                (c.carbon_source == cs)
                and (c.species_names[0] == "Oa")
                and (len(c.species_names) == 1)
                and (c.concentration == conc)
                and (c.signal == "OD")
            ):
                for x, y in zip(c.xs, c.ys):
                    mask = [True if xi < 27 else False for xi in x]
                    x = x[mask]
                    y = y[mask]
                    fig.add_trace(
                        go.Scatter(
                            x=x,
                            y=y,
                            mode="lines",
                            name=f"{cs} {conc} mM",
                            line=dict(color=colors[cs]),
                            marker=dict(size=6),
                            showlegend=False,
                        )
                    )
    for cs in carbon_sources:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines+markers",
                name=f"{cs}",
                line=dict(color=colors[cs]),
                marker=dict(size=6),
            )
        )
    fig.update_layout(
        title=f"Growth Curves of stocks at {', '.join(c_to_mM[conc_init])}",
        xaxis=dict(title="Time [h]", ticks="inside"),
        yaxis=dict(title="OD600", ticks="inside"),
        legend=dict(
            title="Carbon Source",
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(0,0,0,0)",
        ),
    )
    fig = style_plot(fig, line_thickness=2, font_size=12)
    fig.write_image(f"plots/growth_curves_{conc_init}mM.svg")


plot_growth_curves_stocks()
