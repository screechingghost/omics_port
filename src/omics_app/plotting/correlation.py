"""
Port of generate_correlation_plot_advanced (app_12-02.R, line 8148) --
only the default "Clustered Heatmap (with dendrogram)" viz_type is
ported (a hierarchically-reordered sample x sample correlation
heatmap). The other nine of R's 10 (Type-1 pairs panels, shade
corrplot, circle, number, ComplexHeatmap-style, triangular, ellipse
corrgram, network, dendrogram-only) aren't ported -- see
ui/visualization.py's Visualization dropdown, which only offers
"Clustered Heatmap" for now. Each is a genuinely different plot type
(a scatterplot matrix, a network graph, ...), not a styling variant of
this one, so none of them can be gotten "for free" from this function.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform


def build_correlation_figure(
    intensity_matrix: pd.DataFrame,
    method: str = "pearson",
    show_values: bool = True,
) -> go.Figure:
    """intensity_matrix: rows = proteins, columns = samples."""
    if intensity_matrix.shape[1] < 2:
        fig = go.Figure()
        fig.add_annotation(text="Insufficient data for correlation", showarrow=False)
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    cor_matrix = intensity_matrix.corr(method=method)
    samples = cor_matrix.columns.tolist()

    dist = 1 - cor_matrix.to_numpy()
    np.fill_diagonal(dist, 0)
    dist = (dist + dist.T) / 2  # guard against asymmetry from floating-point drift
    dist = np.clip(dist, 0, None)

    order = list(range(len(samples)))
    if len(samples) > 2:
        try:
            z = linkage(squareform(dist, checks=False), method="complete")
            order = list(leaves_list(z))
        except Exception:
            pass

    ordered = cor_matrix.iloc[order, order]
    labels = [samples[i] for i in order]
    values = ordered.to_numpy()

    fig = go.Figure(
        go.Heatmap(
            z=values,
            x=labels,
            y=labels,
            colorscale="RdBu_r",
            zmin=-1,
            zmax=1,
            text=np.round(values, 2) if show_values else None,
            texttemplate="%{text}" if show_values else None,
            colorbar=dict(title=method.capitalize()),
        )
    )
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(
        title=f"Sample Correlation Heatmap ({method.capitalize()})",
        template="plotly_white",
        xaxis_tickangle=-45,
        margin=dict(t=50, b=80),
    )
    return fig
