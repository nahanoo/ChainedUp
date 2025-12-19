from experiment import Species, Experiment
import plotly.graph_objects as go
from style import colors, style_plot


def plot_chemostat_growth_curves():
    e = Experiment(d="../../data/251214_at_chemostat_batch_growth/metaod/")
    e.build_conditions()

    conditions = [
        "At chemostat",
        "Succinate",
        "Glucose",
        "Succinate+Glucose",
        "Succinate+Glucose Outflow",
    ]

    adjust_od = {
        "Succinate": 0.16,
        "Glucose": 0.15,
        "Succinate+Glucose": 0.16,
        "Succinate+Glucose Outflow": 0.55,
    }

    fig = go.Figure()
    for c in e.conditions:
        if (c.species_names[0] in conditions) and (c.carbon_source in conditions):
            for x, y in zip(c.xs, c.ys):
                y0 = adjust_od[c.carbon_source] - y[0]
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y + y0,
                        mode="lines",
                        name=f"{c.carbon_source}",
                        showlegend=False,
                        line=dict(color=colors[c.carbon_source]),
                    )
                )
    for c in conditions[1:]:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines+markers",
                name=f"{c}",
                line=dict(color=colors[c]),
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
        ),
    )
    fig = style_plot(fig, marker_size=6, line_thickness=2.5, font_size=12)
    fig.write_image("plots/chemostat_growth_curve.svg")


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
    conc_init = 5
    concentrations = [conc_init, conc_init, 2 * conc_init]
    for cs, conc in zip(carbon_sources, concentrations):
        for c in e.conditions:
            if (
                (c.carbon_source == cs)
                and (c.species_names[0] == "At")
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
        ),
    )
    fig = style_plot(fig, line_thickness=2, font_size=12)
    fig.write_image(f"plots/growth_curves_{conc_init}mM.svg")
