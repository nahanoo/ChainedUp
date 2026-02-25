from style import colors, style_plot
from chibio_parser import cfu_parser
import plotly.graph_objects as go
from plotly.subplots import make_subplots


df, _ = cfu_parser("260130_at_oa_chemostat_d_03")
reactor_to_substrate = {
    "M0": "Succinate",
    "M1": "Glucose",
    "M2": "Succinate+Glucose",
    "M3": "Succinate+Glucose Outflow",
}
df["average"] = df["average"].replace(0, None)

fig = make_subplots(
    rows=2,
    cols=2,
    subplot_titles=list(reactor_to_substrate.values()),
    shared_xaxes=True,
    shared_yaxes=True,
    vertical_spacing=0.05,
    horizontal_spacing=0.05,
)


for i, (sp_name, sp) in enumerate(df.groupby("species")):
    for j, (r_name, r) in enumerate(sp.groupby("reactor")):
        row = (j // 2) + 1
        col = (j % 2) + 1
        fig.add_trace(
            go.Scatter(
                x=r["sample_time"],
                y=r["average"],
                mode="lines+markers",
                name=f"{sp_name.capitalize()}",
                line=dict(color=colors[sp_name]),
                marker=dict(symbol="circle", size=8, color=colors[sp_name]),
                error_y=dict(type="data", array=r["stdev"], visible=True),
                legendgroup=sp_name,
                showlegend=True if (row == 1 and col == 1) else False,
            ),
            row=row,
            col=col,
        )
fig.for_each_yaxis(
    lambda yaxis: yaxis.update(type="log", range=[5, 10], dtick=1, title_text="CFU/mL")
)
fig.for_each_xaxis(lambda xaxis: xaxis.update(title_text="Time [h]"))
for i in range(1, 3):
    fig.update_yaxes(title=None, row=i, col=2)
    fig.update_xaxes(title=None, row=1, col=i)

fig = style_plot(fig, line_thickness=2, marker_size=6, font_size=12)
fig.write_image("plots/cfu_chemostat_species.svg")
