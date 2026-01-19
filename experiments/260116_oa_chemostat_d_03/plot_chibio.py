import plotly.graph_objects as go
from chibio_parser import fluorescence_paresr, calibration_csv
from style import style_plot, colors
import pandas as pd
import glob
import plotly.express as px

reactor_map = {
    "M0": "Succinate",
    "M1": "Glucose",
    "M2": "Succinate+Glucose",
    "M3": "Succinate+Glucose Outflow",
}


def get_od():
    df = fluorescence_paresr("260119_oa_chemostat_d_03")
    df = df[(df["exp_time"] >= 0.08) & (df["exp_time"] <= 140)]
    df = calibration_csv("../../data/260119_oa_chemostat_d_03/calibration.csv", df)
    return df


def plot_chibio(df):
    reactors = sorted(df["reactor"].unique())
    fig = go.Figure()
    for reactor in reactors:
        od_r = df[df["reactor"] == reactor]
        fig.add_trace(
            go.Scatter(
                x=od_r["exp_time"],
                y=od_r["od_calibrated"],
                mode="lines",
                name=f"{reactor_map[reactor]}",
                line=dict(color=colors[reactor_map[reactor]]),
                marker=dict(size=6),
            )
        )

    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        yaxis=dict(ticks="inside", title="OD600"),
        title="At in chemostats at dilution rate 0.03 h<sup>-1</sup>",
    )
    return fig


def plot_chibio_dilution():
    df = get_od()
    df.to_csv("dataframes/chibio_od.csv", index=False)
    fig = plot_chibio(df)
    fig.update_layout(
        legend=dict(
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(0,0,0,0)",
        ),
    )
    fig.add_vline(
        x=93,
        line=dict(color="gray", dash="dash"),
        annotation=dict(
            text="Stopped dilution",
            showarrow=False,
        ),
    )
    fig = style_plot(fig, marker_size=5)
    fig.write_image("plots/chibio_od_full.svg")

    df = df[(df["exp_time"] <= 93)]
    fig = plot_chibio(df)
    fig = style_plot(fig, marker_size=5)
    fig.update_layout(
        legend=dict(
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(0,0,0,0)",
        ),
    )
    fig.write_image("plots/chibio_od.svg")


def plot_chibio_succ_gluc():
    df = get_od()
    df.to_csv("dataframes/chibio_od.csv", index=False)
    df = df[(df["reactor"] == "M0") | (df["reactor"] == "M1")]
    df = df[(df["exp_time"] <= 93)]
    fig = plot_chibio(df)
    fig.update_layout(
        legend=dict(
            title="Carbon Source",
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(0,0,0,0)",
        ),
    )
    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        # yaxis=dict(ticks="inside", title="CFU/mL", type="log", range=[6, 10]),
        title="Mono-culture Succinate+Glucose",
        legend=dict(title="Carbon Source", x=0.00, y=0.01),
        height=190,
        width=250,
        showlegend=False,
    )

    fig = style_plot(fig, line_thickness=2, marker_size=6, font_size=14)
    fig.write_image("plots/chibio_od_succinate_glucose.svg")


def plot_chain():
    df = get_od()
    df.to_csv("dataframes/chibio_od.csv", index=False)
    df = df[(df["reactor"] == "M2") | (df["reactor"] == "M3")]
    df = df[(df["exp_time"] <= 93)]
    fig = plot_chibio(df)
    fig.update_layout(
        legend=dict(
            title="Carbon Source",
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
        ),
    )
    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        # yaxis=dict(ticks="inside", title="CFU/mL", type="log", range=[6, 10]),
        title="Mono-culture Succinate+Glucose",
        legend=dict(title="Carbon Source", x=0.00, y=0.01),
        height=190,
        width=250,
        showlegend=False,
    )

    fig = style_plot(fig, line_thickness=2, marker_size=5, font_size=14)
    fig.write_image("plots/chibio_od_succinate_glucose_chain.svg")


def parse_plate_reader(dir):
    fs = sorted(
        glob.glob(
            f"{dir}/*absorbance*.xlsx",
            include_hidden=False,
        )
    )
    od = pd.DataFrame(
        columns=[
            "time",
            "Succinate",
            "Succinate_std",
            "Glucose",
            "Glucose_std",
            "Succinate+Glucose",
            "Succinate+Glucose_std",
            "Succinate+Glucose Outflow",
            "Succinate+Glucose Outflow_std",
        ]
    )
    for i, f in enumerate(fs):
        df = pd.read_excel(
            f,
            skiprows=24,
            nrows=1,
            index_col=0,
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        )
        od_s, od_s_std = df[[1, 2, 3]].to_numpy().mean(), df[[1, 2, 3]].to_numpy().std()
        od_g, od_g_std = df[[4, 5, 6]].to_numpy().mean(), df[[4, 5, 6]].to_numpy().std()
        od_sg, od_sg_std = (
            df[[7, 8, 9]].to_numpy().mean(),
            df[[7, 8, 9]].to_numpy().std(),
        )
        od_sgo, od_sgo_std = (
            df[[10, 11, 12]].to_numpy().mean(),
            df[[10, 11, 12]].to_numpy().std(),
        )
        od.loc[len(od)] = [
            i,
            od_s,
            od_s_std,
            od_g,
            od_g_std,
            od_sg,
            od_sg_std,
            od_sgo,
            od_sgo_std,
        ]
    return od


def plot_plate_reader():
    df = parse_plate_reader("../../data/260119_oa_chemostat_d_03/plate_reader")
    df.to_csv("dataframes/plate_reader_od.csv", index=False)
    df.replace(0, None, inplace=True)
    fig = go.Figure()
    for condition in [
        "Succinate",
        "Glucose",
        "Succinate+Glucose",
        "Succinate+Glucose Outflow",
    ]:
        fig.add_trace(
            go.Scatter(
                x=df["time"],
                y=df[condition],
                error_y=dict(
                    type="data",
                    array=df[f"{condition}_std"],
                    visible=True,
                    thickness=1.5,
                    width=3,
                ),
                mode="lines+markers",
                name=condition,
                line=dict(color=colors[condition]),
                marker=dict(size=6),
            )
        )
    fig.update_layout(
        xaxis=dict(ticks="inside", title="Time [h]"),
        yaxis=dict(ticks="inside", title="OD600"),
        legend=dict(
            title="Carbon Source",
            x=0.99,
            y=0.6,
            xanchor="right",
            yanchor="middle",
        ),
        title="At in chemostats at dilution rate 0.03 h<sup>-1</sup> (Plate Reader)",
    )
    fig = style_plot(fig, line_thickness=3, marker_size=5)
    fig.write_image("plots/plate_reader_od.svg")


# plot_chibio_dilution()
# plot_chibio_succ_gluc()
# plot_chibio_succ_gluc()
plot_plate_reader()
