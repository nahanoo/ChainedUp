import pandas as pd
import plotly.graph_objects as go
from style import *

df = pd.read_csv("data.csv")
df["average"] = df["average"].replace(0, None)


species = [["at"], ["oa"], ["at", "oa"]]
reactor_to_substrate = {
    "M0": "Succinate",
    "M1": "Glucose",
    "M2": "Succinate Glucose",
    "M3": "Succinate Glucose Outflow",
}

for exp_name, exp in df.groupby("comment"):
    for r_name, r in exp.groupby("reactor"):
        fig = go.Figure()
        for sp_name, sp in r.groupby("species"):
            fig.add_trace(
                go.Scatter(
                    x=sp["sample_time"],
                    y=sp["average"],
                    mode="lines+markers",
                    name=f"{sp_name.capitalize()}",
                    marker=dict(
                        symbol="circle",
                        size=8,
                        color=colors[sp_name],
                        line=dict(color="black", width=1.2),
                    ),
                    showlegend=False,
                )
            )
        fig.update_layout(
            title=f"{exp_name} {reactor_to_substrate[r_name]}",
            xaxis=dict(title="Time (h)", range=[0, 100], dtick=25),
            yaxis=dict(
                type="log",
                title="CFU/mL",
                exponentformat="power",
                range=[6, 11],
                dtick=2,
            ),
            legend=dict(
                yanchor="top",
                y=0.95,
                xanchor="left",
                x=0.05,
                bgcolor="rgba(255,255,255,0)",
            ),
            width=300,
            height=200,
        )
        fig = style_plot(fig, line_thickness=3, marker_size=10, font_size=12)
        fig.write_image(f"plots/cfu_{exp_name}_{r_name}.svg")
