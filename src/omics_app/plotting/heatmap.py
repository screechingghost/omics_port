"""
Port of generate_heatmap_advanced (app_12-02.R, line 6606) -- 2-group
path only (the ANOVA/multi-group branch, R lines 6700-6841, isn't
ported -- same known gap as everywhere else in this codebase).

Known simplifications vs. R:
- R uses ComplexHeatmap with k-means row/column splitting (row_km /
  column_km) layered on top of the requested hierarchical clustering,
  producing visually separated blocks. This port does plain
  hierarchical reordering only -- no k-means split blocks.
- Distance options: euclidean, manhattan, canberra, maximum, and the
  three correlation-based ones (pearson/spearman/kendall) are
  supported. R's remaining two (binary, minkowski) fall back to
  euclidean here -- they're rarely meaningful for continuous
  log-abundance data anyway.
- "kendall" distance is approximated with a Spearman-style rank
  transform + correlation distance (true all-pairs Kendall's tau is
  O(n^2 log n) per pair and not worth it for a clustering order that's
  already an approximation of R's ComplexHeatmap layout).
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

_METRIC_MAP = {
    "euclidean": "euclidean",
    "manhattan": "cityblock",
    "canberra": "canberra",
    "maximum": "chebyshev",
    "pearson": "correlation",
    "spearman": "correlation",
    "kendall": "correlation",
}
_RANK_BASED = {"spearman", "kendall"}


def _cluster_order(matrix: np.ndarray, distance: str, method: str) -> list[int]:
    """matrix: rows to cluster (observations x features)."""
    if matrix.shape[0] < 3:
        return list(range(matrix.shape[0]))
    if distance in _RANK_BASED:
        matrix = pd.DataFrame(matrix).rank(axis=1).to_numpy()
    metric = _METRIC_MAP.get(distance, "euclidean")
    linkage_method = "average" if method in ("ward.D", "ward.D2") else method
    try:
        dist = pdist(matrix, metric=metric)
        if not np.isfinite(dist).all():
            finite = dist[np.isfinite(dist)]
            dist = np.nan_to_num(dist, nan=(finite.max() if finite.size else 1.0))
        z = linkage(dist, method=linkage_method)
        return list(leaves_list(z))
    except Exception:
        return list(range(matrix.shape[0]))


def select_heatmap_rows(
    df: pd.DataFrame, sig_mask: pd.Series, mode: str, top_n: int
) -> pd.DataFrame:
    """Port of select_heatmap_rows (R lines 6658-6695)."""
    sig = df[sig_mask.reindex(df.index, fill_value=False)]
    if sig.empty or "logFC" not in sig.columns:
        return sig.head(top_n) if mode != "all" else sig

    mode = (mode or "all").lower()
    if mode == "all":
        return sig

    up_ids = sig[sig["logFC"] > 0].sort_values("logFC", ascending=False).head(int(top_n)).index
    down_ids = sig[sig["logFC"] < 0].sort_values("logFC").head(int(top_n)).index
    if mode == "up":
        rows = list(up_ids)
    elif mode == "down":
        rows = list(down_ids)
    elif mode == "both":
        rows = list(up_ids) + [i for i in down_ids if i not in set(up_ids)]
    else:
        rows = list(sig.index)
    return sig.loc[rows] if rows else sig.head(int(top_n))


def build_heatmap_figure(
    df: pd.DataFrame,
    sig_mask: pd.Series,
    sample_columns: list[str],
    sample_groups: dict[str, str],
    heatmap_mode: str = "all",
    top_n: int = 50,
    scale_data: bool = True,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    row_distance: str = "euclidean",
    col_distance: str = "euclidean",
    row_method: str = "complete",
    col_method: str = "complete",
    show_row_names: bool = True,
    color_mapping: dict[str, str] | None = None,
) -> go.Figure:
    selected = select_heatmap_rows(df, sig_mask, heatmap_mode, top_n)
    if selected.empty:
        fig = go.Figure()
        fig.add_annotation(text="No significant proteins to display", showarrow=False)
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    matrix = selected[sample_columns].to_numpy(dtype=float)
    if scale_data:
        mean = np.nanmean(matrix, axis=1, keepdims=True)
        sd = np.nanstd(matrix, axis=1, ddof=1, keepdims=True)
        sd[sd == 0] = 1.0
        matrix = (matrix - mean) / sd
    matrix = np.nan_to_num(matrix)

    row_order = (
        _cluster_order(matrix, row_distance, row_method)
        if cluster_rows and matrix.shape[0] > 2
        else list(range(matrix.shape[0]))
    )
    col_order = (
        _cluster_order(matrix.T, col_distance, col_method)
        if cluster_cols and matrix.shape[1] > 2
        else list(range(matrix.shape[1]))
    )

    matrix = matrix[np.ix_(row_order, col_order)]
    row_labels = (
        [str(selected.index[i]) for i in row_order] if show_row_names else ["" for _ in row_order]
    )
    col_names = [sample_columns[i] for i in col_order]

    fig = go.Figure(
        go.Heatmap(
            z=matrix,
            x=col_names,
            y=row_labels,
            colorscale="RdYlBu_r",
            zmid=0,
            colorbar=dict(title="Scaled<br>Abundance"),
        )
    )
    fig.update_xaxes(tickangle=-45)
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(
        title=f"Heatmap ({len(selected)} proteins, mode={heatmap_mode})",
        template="plotly_white",
        margin=dict(t=50, b=80),
    )
    return fig
