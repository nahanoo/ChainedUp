import pandas as pd
import plotly.graph_objects as go
from style import style_plot, colors
import plotly.express as px
import numpy as np


def plot_yield():
    df = pd.read_csv("growth_yields.csv", dtype={"max_od": float})
    substrates = ["Succinate", "Glucose", "Succinate+Glucose"]

    for substrate in substrates:
        df_sub = df[df["substrate"] == substrate]

        fig = go.Figure()
        for s in ["At", "Oa", "At+Oa"]:
            df_species = df_sub[df_sub["species"] == s]
            df_species.index = range(len(df_species))
            if substrate != "Succinate+Glucose":
                for i in ["max_od_r1", "max_od_r2", "max_od_r3"]:
                    fig.add_trace(
                        go.Scatter(
                            x=df_species["m_conc"].astype(float),
                            y=df_species[i],
                            marker=dict(color=colors[s]),
                            mode="markers",
                            name=s,
                            showlegend=False,
                        )
                    )
                fig.update_xaxes(
                    title=f"Concentration (mM {substrate})",
                    tickmode="array",
                    tickvals=df_species["m_conc"].astype(float),
                    ticktext=df_species["m_conc"].astype(float),
                )
            else:
                df_species.index = range(len(df_species))

                labels = []
                total_conc = []
                for i, conc in enumerate(df_species["m_conc"]):
                    conc = conc.replace("(", "").replace(")", "")
                    m_succ, m_gluc = conc.split(", ")
                    total_conc.append(float(m_succ) + float(m_gluc))
                    labels.append(f"{m_succ} mM, {m_gluc} mM")
                for i in ["max_od_r1", "max_od_r2", "max_od_r3"]:
                    fig.add_trace(
                        go.Scatter(
                            x=total_conc,
                            y=df_species[i],
                            marker=dict(color=colors[s]),
                            mode="markers",
                            name=s,
                            showlegend=False,
                        )
                    )
                fig.update_xaxes(
                    tickmode="array",
                    tickvals=total_conc,
                    ticktext=labels,
                    title="Concentration (mM Succinate, mM Glucose)",
                )
        for s in ["At", "Oa", "At+Oa"]:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    name=s,
                    marker=dict(color=colors[s], size=8),
                )
            )

        fig.update_layout(
            title=f"Growth Yield on {substrate}",
            yaxis=dict(title="Yield (OD600/mM)", ticks="inside"),
            legend=dict(title="Species", x=0.01, y=0.99, yanchor="top"),
        )

        fig = style_plot(fig, font_size=12, marker_size=8)
        fig.write_image(f"plots/yield_{substrate.replace('+', '_').lower()}.svg")


def compare_yields():
    df = pd.read_csv("growth_yields.csv", dtype={"max_od": float})
    df_g = df[df["substrate"] == "Glucose"]
    df_g.index = range(len(df_g))
    df_s = df[df["substrate"] == "Succinate"]
    df_s.index = range(len(df_s))
    df_gs = df[df["substrate"] == "Succinate+Glucose"]
    df_gs.index = range(len(df_gs))

    def viridis_palette_plotly(cs):
        xs = np.linspace(0.12, 0.88, len(cs))
        cols = px.colors.sample_colorscale("Viridis", xs)
        return dict(zip(cs, cols))

    palette = viridis_palette_plotly(df_gs["m_conc"].unique())
    fig = go.Figure()
    dash = {"At": "solid", "Oa": "dash", "At+Oa": "dot"}
    for i in range(len(df_g)):
        od_g, od_s, od_gs = [], [], []
        for j in ["max_od_r1", "max_od_r2", "max_od_r3"]:
            od_g.append(df_g.loc[i, j])
            od_s.append(df_s.loc[i, j])
            od_gs.append(df_gs.loc[i, j])
        x = ["Succinate", "Glucose", "Succinate+Glucose"]
        y = [np.mean(od_s), np.mean(od_g), np.mean(od_gs)]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                marker_color=palette[df_gs.loc[i, "m_conc"]],
                line=dict(dash=dash[df_gs.loc[i, "species"]]),
                showlegend=False,
                mode="lines+markers",
                error_y=dict(
                    type="data", array=[np.std(od_g), np.std(od_s), np.std(od_gs)]
                ),
            )
        )
    for c in df_gs["m_conc"].unique():
        c_glucose, c_succinate = c.replace("(", "").replace(")", "").split(", ")
        label = f"{c_glucose} mM Glucose<br>{c_succinate} mM Succinate"
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                name=label,
                marker=dict(color=palette[c], size=8),
            )
        )
    for s in ["At", "Oa", "At+Oa"]:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines+markers",
                name=s,
                line=dict(dash=dash[s], color="black"),
                marker=dict(color="black", size=8),
            )
        )
    fig.update_layout(
        yaxis=dict(title="Maximum OD600"),
        xaxis=dict(title="Carbon Source"),
    )

    fig = style_plot(fig, font_size=12, line_thickness=2, marker_size=8)
    fig.write_image("plots/yield_comparision.svg")
