"""
Port of generate_boxplot_advanced (app_12-02.R, line 7491) -- 2-group
path only. Of R's 4 plot_style options (raincloud/violin/traditional/
jitter), only "violin" and "traditional" are ported. "raincloud" needs
ggdist's half-eye density geometry, which has no direct Plotly
equivalent worth approximating; "jitter" is a straightforward addition
(traditional boxplot + jittered points) left for later if it turns out
to matter.
"""

import pandas as pd
import plotly.graph_objects as go


def build_boxplot_figure(
    intensity_matrix: pd.DataFrame,
    sample_groups: dict[str, str],
    plot_style: str = "violin",
    color_mapping: dict[str, str] | None = None,
) -> go.Figure:
    """intensity_matrix: rows = proteins, columns = samples (log2,
    filtered, imputed, normalized data)."""
    if intensity_matrix.shape[1] == 0:
        fig = go.Figure()
        fig.add_annotation(text="No sample data available", showarrow=False)
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    long = intensity_matrix.melt(var_name="Sample", value_name="Expression")
    long["Group"] = long["Sample"].map(sample_groups).fillna("Unknown")

    fig = go.Figure()
    seen_groups: set[str] = set()
    for sample in intensity_matrix.columns:
        sub = long[long["Sample"] == sample]
        group = sample_groups.get(sample, "Unknown")
        color = (color_mapping or {}).get(group)
        show_legend = group not in seen_groups
        seen_groups.add(group)

        if plot_style == "traditional":
            fig.add_trace(
                go.Box(
                    y=sub["Expression"],
                    name=sample,
                    marker_color=color,
                    legendgroup=group,
                    showlegend=show_legend if show_legend else False,
                )
            )
        else:  # violin (default)
            fig.add_trace(
                go.Violin(
                    y=sub["Expression"],
                    name=sample,
                    line_color=color,
                    legendgroup=group,
                    showlegend=False,
                    box_visible=True,
                    meanline_visible=True,
                )
            )
        if show_legend:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(color=color, size=10),
                    name=group,
                    legendgroup=group,
                    showlegend=True,
                )
            )

    fig.update_layout(
        title="Intensity Distribution by Sample",
        yaxis_title="Log2 Abundance",
        xaxis_title="Samples",
        template="plotly_white",
        xaxis_tickangle=-45,
        margin=dict(t=50, b=80),
    )
    return fig
