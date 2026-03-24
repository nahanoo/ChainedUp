import plotly.graph_objects as go
import pandas as pd
import numpy as np
from chibio_parser import *
from style import *
import conditions as cond
from plotly.express.colors import qualitative


def add_lod(fig, min_y, y_lod=1e-3):
    if min_y < 1e-2:
        fig.add_hline(
            y=y_lod,
            line=dict(color=qualitative.Bold[10], dash="dash", width=2),
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


def make_resource_figure(df, wells, xs, ys, title, x_range, dtick):
    selected = [
        (well, x) for well, x in zip(wells, xs) if x_range[0] <= x <= x_range[1]
    ]
    wells_sel = [well for well, _ in selected]
    xs_sel = [x for _, x in selected]

    fig = go.Figure()
    min_ys = []

    for ycol in ys:
        key = ycol.split("_")[0].capitalize()
        yvals = [df.loc[well, ycol] for well in wells_sel]
        min_ys.append(min(yvals))

        fig.add_trace(
            go.Scatter(
                x=xs_sel,
                y=yvals,
                mode="markers+lines",
                name=key[:4],
                marker=dict(
                    color=colors[key],
                    line=dict(color="black", width=1.2),
                ),
            )
        )

    fig.update_layout(
        width=200,
        height=180,
        legend=dict(
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=0.95,
            bgcolor="rgba(255,255,255,0)",
        ),
        xaxis=dict(title="Time (h)", range=x_range, dtick=dtick),
        yaxis=dict(
            title="Area",
            type="log",
            exponentformat="power",
        ),
        title=title,
        showlegend=False,
    )

    fig = style_plot(fig, marker_size=8, line_thickness=3, font_size=11)
    add_lod(fig, min(min_ys))
    return fig


def plot_resources():
    df = pd.read_csv("normalized_data.csv", index_col="id")

    conds = [
        cond.OA_SUCCINATE,
        cond.OA_GLUCOSE,
        cond.OA_SUCCINATE_GLUCOSE,
    ]

    titles = [
        "Oa Succinate",
        "Oa Glucose",
        "Oa Succinate Glucose",
    ]

    xs_list = [
        cond.OA_TIMEPOINTS,
        cond.OA_TIMEPOINTS,
        cond.OA_TIMEPOINTS,
    ]

    ys_list = [
        ["succinate_area_normalized"],
        ["glucose_area_normalized"],
        ["succinate_area_normalized", "glucose_area_normalized"],
    ]

    fnames = [
        "oa_succinate.svg",
        "oa_glucose.svg",
        "oa_succinate_glucose.svg",
    ]

    zoom_specs = {
        "oa_succinate_glucose.svg": {
            "zoom_fname": "oa_succinate_glucose_0_60.svg",
            "x_range": [0, 60],
            "dtick": 20,
        }
    }

    for wells, title, xs, ys, fname in zip(conds, titles, xs_list, ys_list, fnames):
        fig = make_resource_figure(
            df=df,
            wells=wells,
            xs=xs,
            ys=ys,
            title=title,
            x_range=[0, 100],
            dtick=25,
        )
        fig.write_image(f"plots/resources/{fname}")

        if fname in zoom_specs:
            spec = zoom_specs[fname]
            fig_zoom = make_resource_figure(
                df=df,
                wells=wells,
                xs=xs,
                ys=ys,
                title=f"{title} (0–60 h)",
                x_range=spec["x_range"],
                dtick=spec["dtick"],
            )
            fig_zoom.write_image(f"plots/resources/{spec['zoom_fname']}")


def plot_cfus():
    df = pd.read_csv("cfus.csv")
    df = df[df["species"] == "oa"]
    df = df[df["reactor"].isin(["Succinate", "Glucose", "Succinate+Glucose"])]

    for cs_name, cs in df.groupby("reactor"):
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=cs["sample_time"],
                y=cs["average"],
                mode="lines+markers",
                name=cs_name,
                marker=dict(
                    color=colors[cs_name],
                    line=dict(color="black", width=1.2),
                ),
            )
        )

        fig.update_layout(
            xaxis=dict(title="Time (h)", range=[0, 100], dtick=25),
            yaxis=dict(
                title="CFU/mL",
                type="log",
                exponentformat="power",
                range=[6, 10],
            ),
            title=f"CFUs/mL Oa in {cs_name.lower()}",
            width=200,
            height=180,
            showlegend=False,
        )

        fig = style_plot(fig, line_thickness=2, marker_size=8, font_size=11)
        fname = f"cfu_{cs_name.lower()}.svg"
        fig.write_image(f"plots/cfus/{fname}")


def plot_od_succinate_glucose():
    df = fluorescence_paresr("260119_oa_chemostat_d_03")
    df = df[(df["exp_time"] >= 0.08) & (df["exp_time"] <= 140)]
    df = calibration_csv("../../data/260119_oa_chemostat_d_03/calibration.csv", df)
    df = df[df["reactor"] == "M2"]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df["exp_time"][::20],
            y=df["od_calibrated"][::20],
            mode="lines",
            line=dict(color=colors["Succinate+Glucose"]),
        )
    )
    fig.update_layout(
        xaxis=dict(title="Time (h)"),
        yaxis=dict(title="OD<sub>600</sub>"),
        title="OD<sub>600</sub> in succinate+glucose",
        width=200,
        height=180,
        showlegend=False,
    )
    fig = style_plot(fig, line_thickness=3, marker_size=8, font_size=11)
    fig.write_image("plots/cfus/od_succinate_glucose.svg")


plot_cfus()
