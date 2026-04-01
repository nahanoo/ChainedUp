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
    bottom_margin=30,
    font_size=11,
    line_thickness=3,
):
    """Apply consistent styling to Plotly figures."""
    fig.update_layout(
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#FFFFFF",
        font=dict(size=font_size, color="black"),
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=bottom_margin),
        hoverlabel=dict(font_size=font_size),
    )

    for trace in fig.data:
        trace_type = getattr(trace, "type", None)

        # Marker size only for traces that support it
        if trace_type in {"scatter", "scattergl"}:
            if hasattr(trace, "marker") and trace.marker is not None:
                trace.marker.size = marker_size

        # Line width for line-based traces
        if trace_type in {"scatter", "scattergl"}:
            if hasattr(trace, "line") and trace.line is not None:
                trace.line.width = line_thickness

        # Error bar thickness if present
        if hasattr(trace, "error_y") and trace.error_y is not None:
            try:
                trace.error_y.thickness = line_thickness
            except Exception:
                pass

    # Style annotations if present
    if fig.layout.annotations:
        for ann in fig.layout.annotations:
            if ann.font is None:
                ann.font = {}
            ann.font.size = font_size
            ann.font.color = "black"

    # Style title if present
    if fig.layout.title is not None:
        if fig.layout.title.font is None:
            fig.layout.title.font = {}
        fig.layout.title.font.size = font_size
        fig.layout.title.font.color = "black"

    # Style legend title if present
    if fig.layout.legend is not None and fig.layout.legend.title is not None:
        if fig.layout.legend.title.font is None:
            fig.layout.legend.title.font = {}
        fig.layout.legend.title.font.size = font_size
        fig.layout.legend.title.font.color = "black"

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
        ticks="inside",
    )

    fig.update_xaxes(
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
        ticks="inside",
    )

    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )

    return fig
