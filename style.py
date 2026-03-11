from plotly.express.colors import qualitative as qualitative

width, height = 180, 180
colors = {
    "At": qualitative.Bold[2],
    "at": qualitative.Bold[2],
    "Oa": qualitative.Bold[4],
    "oa": qualitative.Bold[4],
    "succinate": qualitative.Bold[5],
    "Succinate": qualitative.Bold[5],
    "glucose": qualitative.Bold[6],
    "Glucose": qualitative.Bold[6],
    "succinate glucose": qualitative.Bold[0],
    "Succinate Glucose": qualitative.Bold[0],
    "Succinate+Glucose": qualitative.Bold[3],
    "Ribose": qualitative.Bold[2],
}


def style_plot(
    fig,
    marker_size=3,
    top_margin=30,
    left_margin=30,
    right_margin=30,
    buttom_margin=30,
    font_size=11,
    line_thickness=3,
):
    """Style function for figures setting fot size and true black color."""
    fig.update_layout(
        {
            "plot_bgcolor": "#FFFFFF",
            "paper_bgcolor": "#FFFFFF",
        },
        font={"size": font_size, "color": "black"},
    )
    for d in fig["data"]:
        try:
            d["marker"]["size"] = marker_size
        except KeyError:
            pass
        try:
            d["line"]["width"] = line_thickness
        except KeyError:
            pass
        try:
            d["error_y"]["thickness"] = line_thickness
        except KeyError:
            pass
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = font_size
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = font_size
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = font_size
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"

    fig.update_layout(
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=buttom_margin),
        hoverlabel=dict(font_size=font_size),
    )
    gridline_width = 0.2
    fig.update_yaxes(
        title_standoff=0,
        gridcolor="gray",
        zeroline=False,
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=0.5,
        showline=True,
        mirror=True,
        linecolor="black",
        linewidth=0.5,
        tickcolor="black",
        tickwidth=0.5,
    )
    fig.update_xaxes(
        title_standoff=0,
        gridcolor="gray",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=0.5,
        showline=True,
        mirror=True,
        linecolor="black",
        linewidth=0.5,
        zeroline=False,
        tickcolor="black",
        tickwidth=0.5,
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_xaxis(lambda axis: axis.update(ticks="inside"))
    fig.for_each_yaxis(lambda axis: axis.update(ticks="inside"))
    return fig
