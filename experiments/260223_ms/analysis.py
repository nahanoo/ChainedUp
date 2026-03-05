import conditions as cond
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from style import *
from plotly.express.colors import qualitative as qualitative

data = pd.read_csv("normalized_data.csv", index_col="id")


def plot_calibration_curves():
    cal_p1_succinate = np.array(
        [data.loc[well, "succinate_area_normalized"] for well in cond.CAL_P1_SUCCINATE]
    )
    cal_p1_glucose = np.array(
        [data.loc[well, "glucose_area_normalized"] for well in cond.CAL_P1_GLUCOSE]
    )
    cal_p2_succinate = np.array(
        [data.loc[well, "succinate_area_normalized"] for well in cond.CAL_P2_SUCCINATE]
    )
    cal_p2_glucose = np.array(
        [data.loc[well, "glucose_area_normalized"] for well in cond.CAL_P2_GLUCOSE]
    )
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=cond.SUCCINATE_CONCENTRATIONS,
            y=cal_p1_succinate,
            mode="markers",
            name="Succ",
            marker=dict(
                color=colors["Succinate"],
                line=dict(color="black", width=1.2),
            ),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=cond.GLUCOSE_CONCENTRATIONS,
            y=cal_p1_glucose,
            mode="markers",
            name="Gluc",
            marker=dict(color=colors["Glucose"], line=dict(color="black", width=1.2)),
        )
    )
    fig.update_layout(
        xaxis=dict(
            type="log",
            exponentformat="power",
            title="Concentration (mM)",
            range=[-4, 1],
            dtick=1,
        ),
        yaxis=dict(
            exponentformat="power", type="log", title="Area", range=[-5, 2], dtick=2
        ),
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="left",
            x=0.05,
            bgcolor="rgba(255,255,255,1)",
        ),
        width=300,
        height=200,
        title="Calibration Curves P1",
    )
    fig = style_plot(fig, marker_size=10, line_thickness=3, font_size=14)
    fig.write_image("plots/calibration_curves_p1.svg")

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=cond.SUCCINATE_CONCENTRATIONS,
            y=cal_p2_succinate,
            mode="markers",
            name="Succ",
            marker=dict(
                color=colors["Succinate"],
                line=dict(color="black", width=1.2),
            ),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=cond.GLUCOSE_CONCENTRATIONS,
            y=cal_p2_glucose,
            mode="markers",
            name="Gluc",
            marker=dict(color=colors["Glucose"], line=dict(color="black", width=1.2)),
        )
    )

    fig.update_layout(
        xaxis=dict(
            type="log",
            exponentformat="power",
            title="Concentration (mM)",
            range=[-4, 1],
            dtick=1,
        ),
        yaxis=dict(
            exponentformat="power", type="log", title="Area", range=[-5, 2], dtick=2
        ),
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="left",
            x=0.05,
            bgcolor="rgba(255,255,255,1)",
        ),
        width=300,
        height=200,
        title="Calibration Curves P2",
    )
    fig = style_plot(fig, marker_size=10, line_thickness=3, font_size=14)
    fig.write_image("plots/calibration_curves_p2.svg")


def plot_all_conditions(cond):
    conds = [
        cond.AT_SUCCINATE,
        cond.AT_GLUCOSE,
        cond.OA_SUCCINATE,
        cond.OA_GLUCOSE,
        cond.ATOA_SUCCINATE,
        cond.ATOA_GLUCOSE,
        cond.AT_SUCCINATE_GLUCOSE,
        cond.ATC2_SUCCINATE_GLUCOSE,
        cond.OA_SUCCINATE_GLUCOSE,
        cond.OAC2_SUCCINATE_GLUCOSE,
        cond.ATOA_SUCCINATE_GLUCOSE,
        cond.ATOAC2_SUCCINATE_GLUCOSE,
    ]

    titles = [
        "At Succinate (Glucose)",
        "At Glucose (2x Glucose)",
        "Oa Succinate",
        "Oa Glucose",
        "AtOA Succinate",
        "AtOA Glucose",
        "At Glucose",
        "AtC2 Glucose",
        "Oa Succinate Glucose",
        "OaC2 Succinate Glucose",
        "AtOA Succinate Glucose",
        "AtOA C2 Succinate Glucose",
    ]

    xs_list = [
        cond.AT_TIMEPOINTS,
        cond.AT_TIMEPOINTS,
        cond.OA_TIMEPOINTS,
        cond.OA_TIMEPOINTS,
        cond.AT_OA_TIMEPOINTS,
        cond.AT_OA_TIMEPOINTS,
        cond.AT_TIMEPOINTS,
        cond.AT_TIMEPOINTS,  # <- likely fix (was AT_OA_TIMEPOINTS)
        cond.OA_TIMEPOINTS,
        cond.OA_TIMEPOINTS,
        cond.AT_OA_TIMEPOINTS,
        cond.AT_OA_TIMEPOINTS,
    ]

    ys_list = [
        ["glucose_area_normalized"],
        ["glucose_area_normalized"],
        ["succinate_area_normalized"],
        ["glucose_area_normalized"],
        ["succinate_area_normalized"],
        ["glucose_area_normalized"],
        ["glucose_area_normalized"],
        ["glucose_area_normalized"],
        ["succinate_area_normalized", "glucose_area_normalized"],
        ["succinate_area_normalized", "glucose_area_normalized"],
        ["succinate_area_normalized", "glucose_area_normalized"],
        ["succinate_area_normalized", "glucose_area_normalized"],
    ]

    fnames = [
        "at_succinate.svg",
        "at_glucose.svg",
        "oa_succinate.svg",
        "oa_glucose.svg",
        "atoa_succinate.svg",
        "atoa_glucose.svg",
        "at_succinate_glucose.svg",
        "atc2_succinate_glucose.svg",
        "oa_succinate_glucose.svg",
        "oac2_succinate_glucose.svg",
        "atoa_succinate_glucose.svg",
        "atoac2_succinate_glucose.svg",
    ]

    for wells, title, xs, ys, fname in zip(conds, titles, xs_list, ys_list, fnames):
        fig = go.Figure()

        for ycol in ys:
            key = ycol.split("_")[0].capitalize()
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=[data.loc[well, ycol] for well in wells],
                    mode="markers+lines",
                    name=key[:4],
                    marker=dict(
                        color=colors[key],
                        line=dict(color="black", width=1.2),
                    ),
                )
            )

        fig.update_layout(
            width=300,
            height=200,
            legend=dict(
                yanchor="top",
                y=0.95,
                xanchor="right",
                x=0.95,
                bgcolor="rgba(255,255,255,0)",
            ),
            xaxis=dict(title="Time (h)", range=[0, 100], dtick=25),
            yaxis=dict(
                title="Area", type="log", exponentformat="power", range=[-5, 2], dtick=2
            ),
            title=title,
            showlegend=False,
        )

        fig = style_plot(fig, marker_size=10, line_thickness=3, font_size=14)
        min_ys = []
        for ycol in ys:
            min_ys.append(min([data.loc[well, ycol] for well in wells]))
        y_lod = 1e-3
        if min(min_ys) < 1e-2:
            fig.add_hline(
                y=y_lod, line=dict(color=qualitative.Bold[10], dash="dash", width=2)
            )
            fig.add_annotation(
                x=10,
                y=np.log10(y_lod),
                text="LOD",
                yanchor="bottom",
                showarrow=False,
                font=dict(size=10, color=qualitative.Bold[10]),
                align="left",
            )
        fig.write_image(f"plots/{fname}")


plot_all_conditions(cond)
plot_calibration_curves()
