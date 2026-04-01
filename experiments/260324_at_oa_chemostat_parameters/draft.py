import os

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

try:
    from style import colors, style_plot
except ImportError:
    from plotly.express.colors import qualitative

    colors = {
        "At": qualitative.Bold[2],
        "Oa": qualitative.Bold[4],
        "Succinate": qualitative.Bold[5],
        "Glucose": qualitative.Bold[6],
        "Succinate+Glucose": qualitative.Bold[3],
    }

    def style_plot(fig, **kwargs):
        return fig


OUTPUT_DIR = "plots/conceptual_allocation"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def allocation_overlap(a_at, a_oa):
    """
    Simple parameter-free proxy:
    overlap is maximal when the two allocations are the same
    and decreases with allocation distance.
    """
    return 1 - np.abs(a_at - a_oa)


def allocation_partitioning(a_at, a_oa):
    return np.abs(a_at - a_oa)


def plot_trait_space():
    a = np.linspace(0, 1, 101)
    A_AT, A_OA = np.meshgrid(a, a)

    overlap = allocation_overlap(A_AT, A_OA)
    partition = allocation_partitioning(A_AT, A_OA)

    fig = go.Figure()

    fig.add_trace(
        go.Heatmap(
            x=a,
            y=a,
            z=overlap,
            colorbar=dict(title="Overlap proxy"),
            colorscale="Viridis",
            showscale=True,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            name="Equal allocation",
            line=dict(color="white", width=2, dash="dash"),
        )
    )

    fig.add_annotation(
        x=0.5,
        y=0.55,
        text="high overlap",
        showarrow=False,
        font=dict(size=13, color="white"),
    )

    fig.add_annotation(
        x=0.1,
        y=0.9,
        text="strong partitioning",
        showarrow=False,
        font=dict(size=13, color="white"),
    )

    fig.add_annotation(
        x=0.9,
        y=0.1,
        text="strong partitioning",
        showarrow=False,
        font=dict(size=13, color="white"),
    )

    fig.add_annotation(
        x=0.05,
        y=-0.08,
        text="At glucose specialist",
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(size=12),
    )
    fig.add_annotation(
        x=0.95,
        y=-0.08,
        text="At succinate specialist",
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(size=12),
    )
    fig.add_annotation(
        x=-0.1,
        y=0.05,
        text="Oa glucose specialist",
        showarrow=False,
        textangle=-90,
        xref="paper",
        yref="paper",
        font=dict(size=12),
    )
    fig.add_annotation(
        x=-0.1,
        y=0.95,
        text="Oa succinate specialist",
        showarrow=False,
        textangle=-90,
        xref="paper",
        yref="paper",
        font=dict(size=12),
    )

    fig.update_layout(
        width=450,
        height=400,
        title="Allocation trait space",
        xaxis=dict(title="At allocation to succinate, a_at", range=[0, 1]),
        yaxis=dict(title="Oa allocation to succinate, a_oa", range=[0, 1]),
    )

    fig = style_plot(fig, font_size=12)
    fig.write_image(os.path.join(OUTPUT_DIR, "trait_space_overlap.svg"))


def plot_overlap_curve():
    delta_a = np.linspace(0, 1, 200)
    overlap = 1 - delta_a
    partition = delta_a

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=delta_a,
            y=overlap,
            mode="lines",
            name="Overlap proxy",
            line=dict(color=colors["Succinate+Glucose"], width=3),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=delta_a,
            y=partition,
            mode="lines",
            name="Partitioning proxy",
            line=dict(color=colors["Oa"], width=3, dash="dash"),
        )
    )

    fig.add_annotation(
        x=0.1,
        y=0.9,
        text="similar allocations → high overlap",
        showarrow=False,
        font=dict(size=12),
    )
    fig.add_annotation(
        x=0.78,
        y=0.82,
        text="divergent allocations → high partitioning",
        showarrow=False,
        font=dict(size=12),
    )

    fig.update_layout(
        width=450,
        height=350,
        title="Allocation distance as a niche-overlap proxy",
        xaxis=dict(title="Allocation distance, |a_at - a_oa|", range=[0, 1]),
        yaxis=dict(title="Proxy value", range=[0, 1]),
        legend=dict(
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255,255,255,0)",
        ),
    )

    fig = style_plot(fig, font_size=12)
    fig.write_image(os.path.join(OUTPUT_DIR, "allocation_distance_curve.svg"))


def plot_exemplar_pairs():
    exemplars = [
        ("both glucose-focused", 0.0, 0.0),
        ("both generalists", 0.5, 0.5),
        ("both succinate-focused", 1.0, 1.0),
        ("partitioned", 0.0, 1.0),
        ("reverse partitioned", 1.0, 0.0),
    ]

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=[x for _, x, _ in exemplars],
            y=[y for _, _, y in exemplars],
            mode="markers+text",
            text=[label for label, _, _ in exemplars],
            textposition="top center",
            marker=dict(
                size=12,
                color=colors["Succinate+Glucose"],
                line=dict(color="black", width=1.2),
            ),
            showlegend=False,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            line=dict(color="gray", width=2, dash="dash"),
            showlegend=False,
        )
    )

    fig.update_layout(
        width=450,
        height=400,
        title="Exemplar strategy pairs",
        xaxis=dict(title="At allocation to succinate, a_at", range=[-0.05, 1.05]),
        yaxis=dict(title="Oa allocation to succinate, a_oa", range=[-0.05, 1.05]),
    )

    fig = style_plot(fig, font_size=12)
    fig.write_image(os.path.join(OUTPUT_DIR, "allocation_exemplar_pairs.svg"))


def plot_allocation_bars():
    pairs = [
        ("both glucose", 0.0, 0.0),
        ("both generalists", 0.5, 0.5),
        ("both succinate", 1.0, 1.0),
        ("partitioned", 0.0, 1.0),
        ("reverse partitioned", 1.0, 0.0),
    ]

    labels = []
    succ_values = []
    gluc_values = []
    bar_colors = []

    for label, a_at, a_oa in pairs:
        labels.extend([f"{label}: At", f"{label}: Oa"])
        succ_values.extend([a_at, a_oa])
        gluc_values.extend([1 - a_at, 1 - a_oa])
        bar_colors.extend([colors["At"], colors["Oa"]])

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            y=labels,
            x=gluc_values,
            name="Glucose allocation",
            orientation="h",
            marker=dict(color=colors["Glucose"]),
        )
    )

    fig.add_trace(
        go.Bar(
            y=labels,
            x=succ_values,
            name="Succinate allocation",
            orientation="h",
            marker=dict(color=colors["Succinate"]),
        )
    )

    fig.update_layout(
        barmode="stack",
        width=700,
        height=450,
        title="Resource-allocation profiles for exemplar strategy pairs",
        xaxis=dict(title="Allocation fraction", range=[0, 1]),
        yaxis=dict(autorange="reversed"),
        legend=dict(
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255,255,255,0)",
        ),
    )

    fig = style_plot(fig, font_size=12)
    fig.write_image(os.path.join(OUTPUT_DIR, "allocation_exemplar_bars.svg"))


def plot_combined_panel():
    a = np.linspace(0, 1, 101)
    A_AT, A_OA = np.meshgrid(a, a)
    overlap = allocation_overlap(A_AT, A_OA)

    delta_a = np.linspace(0, 1, 200)
    overlap_curve = 1 - delta_a
    partition_curve = delta_a

    exemplars = [
        ("both glucose", 0.0, 0.0),
        ("generalists", 0.5, 0.5),
        ("both succinate", 1.0, 1.0),
        ("partitioned", 0.0, 1.0),
        ("rev. partitioned", 1.0, 0.0),
    ]

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "A. Trait space",
            "B. Allocation distance",
            "C. Exemplar pairs",
            "D. Allocation profiles",
        ),
        specs=[
            [{"type": "heatmap"}, {"type": "xy"}],
            [{"type": "xy"}, {"type": "xy"}],
        ],
        vertical_spacing=0.12,
        horizontal_spacing=0.12,
    )

    fig.add_trace(
        go.Heatmap(
            x=a,
            y=a,
            z=overlap,
            colorscale="Viridis",
            showscale=False,
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            line=dict(color="white", width=2, dash="dash"),
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=delta_a,
            y=overlap_curve,
            mode="lines",
            name="Overlap proxy",
            line=dict(color=colors["Succinate+Glucose"], width=3),
            showlegend=True,
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=delta_a,
            y=partition_curve,
            mode="lines",
            name="Partitioning proxy",
            line=dict(color=colors["Oa"], width=3, dash="dash"),
            showlegend=True,
        ),
        row=1,
        col=2,
    )

    fig.add_trace(
        go.Scatter(
            x=[x for _, x, _ in exemplars],
            y=[y for _, _, y in exemplars],
            mode="markers+text",
            text=[label for label, _, _ in exemplars],
            textposition="top center",
            marker=dict(
                size=10,
                color=colors["Succinate+Glucose"],
                line=dict(color="black", width=1.2),
            ),
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode="lines",
            line=dict(color="gray", width=2, dash="dash"),
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    pair_labels = [
        "both glucose",
        "generalists",
        "both succinate",
        "partitioned",
        "rev. partitioned",
    ]
    pair_a_at = [0.0, 0.5, 1.0, 0.0, 1.0]
    pair_a_oa = [0.0, 0.5, 1.0, 1.0, 0.0]

    fig.add_trace(
        go.Bar(
            x=[1 - x for x in pair_a_at],
            y=[f"{lab}: At" for lab in pair_labels],
            orientation="h",
            marker=dict(color=colors["Glucose"]),
            name="Glucose allocation",
            showlegend=False,
        ),
        row=2,
        col=2,
    )
    fig.add_trace(
        go.Bar(
            x=pair_a_at,
            y=[f"{lab}: At" for lab in pair_labels],
            orientation="h",
            marker=dict(color=colors["Succinate"]),
            name="Succinate allocation",
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    fig.add_trace(
        go.Bar(
            x=[1 - x for x in pair_a_oa],
            y=[f"{lab}: Oa" for lab in pair_labels],
            orientation="h",
            marker=dict(color=colors["Glucose"]),
            showlegend=False,
        ),
        row=2,
        col=2,
    )
    fig.add_trace(
        go.Bar(
            x=pair_a_oa,
            y=[f"{lab}: Oa" for lab in pair_labels],
            orientation="h",
            marker=dict(color=colors["Succinate"]),
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    fig.update_xaxes(title_text="a_at", row=1, col=1, range=[0, 1])
    fig.update_yaxes(title_text="a_oa", row=1, col=1, range=[0, 1])

    fig.update_xaxes(title_text="|a_at - a_oa|", row=1, col=2, range=[0, 1])
    fig.update_yaxes(title_text="Proxy value", row=1, col=2, range=[0, 1])

    fig.update_xaxes(title_text="a_at", row=2, col=1, range=[-0.05, 1.05])
    fig.update_yaxes(title_text="a_oa", row=2, col=1, range=[-0.05, 1.05])

    fig.update_xaxes(title_text="Allocation fraction", row=2, col=2, range=[0, 1])
    fig.update_yaxes(autorange="reversed", row=2, col=2)

    fig.update_layout(
        width=1000,
        height=800,
        title="Conceptual allocation plots",
        barmode="stack",
        legend=dict(
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255,255,255,0)",
        ),
    )

    fig = style_plot(fig, font_size=12)
    fig.write_image(os.path.join(OUTPUT_DIR, "allocation_conceptual_panel.svg"))


def main():
    plot_trait_space()
    plot_overlap_curve()
    plot_exemplar_pairs()
    plot_allocation_bars()
    plot_combined_panel()
    print(f"Saved conceptual plots to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
