"""
Port of the enrichment plot renderer (app_12-02.R, generate_enrichment_plot
and friends, ~line 10250 onward) -- only 2 of R's 9 plot_type options
are wired here: Grouped Barplot and Dotplot, matching the same
"default options only" approach used throughout the rest of
plotting/*.py. The other 7 (Bubble Plot, Up vs Down Comparison,
Gene-Pathway Network / cnetplot, Enrichment Map / emapplot, Chord: Up,
Chord: Down) aren't ported -- see ui/enrichment.py's plot-type
dropdown, which shows a placeholder message for the unported options.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

UP_COLOR = "#BB0A1E"
DOWN_COLOR = "royalblue"


def _parse_overlap(overlap: str) -> tuple[int, int]:
    try:
        k, n = str(overlap).split("/")
        return int(k), int(n)
    except Exception:
        return 0, 1


def build_enrichment_barplot(df: pd.DataFrame, top_n: int = 15) -> go.Figure:
    if df is None or df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No enrichment results to display", showarrow=False)
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    data = df.sort_values("PValue").head(top_n).copy()
    data["neglog10p"] = -np.log10(data["PValue"].clip(lower=1e-300))
    data = data.iloc[::-1]  # smallest p-value ends up drawn at the top

    fig = go.Figure()
    for direction, color in [("Up", UP_COLOR), ("Down", DOWN_COLOR)]:
        subset = data[data["Direction"] == direction]
        if subset.empty:
            continue
        fig.add_trace(
            go.Bar(
                x=subset["neglog10p"],
                y=subset["Term"],
                orientation="h",
                name=direction,
                marker_color=color,
                hovertemplate="%{y}<br>-log10(p)=%{x:.2f}<extra></extra>",
            )
        )

    fig.update_layout(
        title=f"Top {min(top_n, len(df))} Enriched Terms",
        xaxis_title="-Log10(p-value)",
        barmode="group",
        template="plotly_white",
        margin=dict(l=300, t=50),
    )
    return fig


def build_enrichment_dotplot(df: pd.DataFrame, top_n: int = 15) -> go.Figure:
    if df is None or df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No enrichment results to display", showarrow=False)
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    data = df.sort_values("PValue").head(top_n).copy()
    parsed = data["Overlap"].apply(_parse_overlap)
    data["k"] = [p[0] for p in parsed]
    data["n"] = [p[1] for p in parsed]
    data["GeneRatio"] = data["k"] / data["n"].replace(0, 1)
    data = data.iloc[::-1]

    max_k = max(int(data["k"].max()), 1)
    fig = go.Figure(
        go.Scatter(
            x=data["GeneRatio"],
            y=data["Term"],
            mode="markers",
            marker=dict(
                size=data["k"],
                sizemode="area",
                sizeref=2.0 * max_k / (40.0**2),
                sizemin=4,
                color=data["PValue"],
                colorscale="Viridis_r",
                showscale=True,
                colorbar=dict(title="p-value"),
            ),
            text=data["Direction"],
            hovertemplate="%{y}<br>GeneRatio=%{x:.2f}<br>%{text}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Top {min(top_n, len(df))} Enriched Terms",
        xaxis_title="Gene Ratio",
        template="plotly_white",
        margin=dict(l=300, t=50),
    )
    return fig
