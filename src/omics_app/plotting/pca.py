"""
Port of generate_pca_plot_advanced (app_12-02.R, line 7049) --
"Biplot (Samples)" only, the default of R's 5 PCA plot types. The
other four -- generate_pca_variable_plot_custom, _scree_plot,
_loading_plot, _cumulative_variance_plot (R lines 7207-7490) -- aren't
ported; see ui/visualization.py, whose "Plot Type" dropdown only
offers Biplot for now.

Known simplification: R's FactoMineR::PCA(scale.unit = TRUE) is
matched here with sklearn's StandardScaler + PCA -- the same
operation (mean-center, unit-variance scale, then SVD) -- but exact PC
sign/orientation can flip between implementations. This is a
well-known, harmless PCA ambiguity (eigenvectors are only defined up
to sign), not a correctness bug; cluster shapes and relative distances
are unaffected.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# 95% confidence ellipse: sqrt of the chi-square(2 df) 0.95 quantile.
_CHI2_95_2DF = 2.4477


def _confidence_ellipse(points: np.ndarray, n_std: float = _CHI2_95_2DF):
    """Returns (x, y) arrays tracing a 95% confidence ellipse for a 2D
    group of points, or None if there are too few points to fit one."""
    if points.shape[0] < 3:
        return None
    cov = np.cov(points, rowvar=False)
    mean = points.mean(axis=0)
    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = np.clip(eigvals, 0, None)
    theta = np.linspace(0, 2 * np.pi, 100)
    circle = np.column_stack([np.cos(theta), np.sin(theta)])
    ellipse = circle @ np.diag(np.sqrt(eigvals) * n_std) @ eigvecs.T + mean
    return ellipse[:, 0], ellipse[:, 1]


def build_pca_figure(
    intensity_matrix: pd.DataFrame,
    sample_groups: dict[str, str],
    pc1: int = 1,
    pc2: int = 2,
    show_ellipses: bool = True,
    show_labels: bool = True,
    color_mapping: dict[str, str] | None = None,
) -> go.Figure:
    """intensity_matrix: rows = proteins, columns = samples (log2,
    filtered, imputed, normalized -- matches R's intensity_matrix)."""
    complete = intensity_matrix.dropna(axis=0)
    if complete.shape[0] < 10 or complete.shape[1] < 3:
        fig = go.Figure()
        fig.add_annotation(text="Insufficient data for PCA", showarrow=False)
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    samples = complete.columns.tolist()
    X = StandardScaler().fit_transform(complete.T.to_numpy())
    max_pc = max(1, min(10, X.shape[0] - 1, X.shape[1]))

    if pc1 > max_pc or pc2 > max_pc:
        fig = go.Figure()
        fig.add_annotation(text=f"Only {max_pc} PCs available", showarrow=False)
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    pca = PCA(n_components=max_pc)
    coords = pca.fit_transform(X)
    var_exp = pca.explained_variance_ratio_ * 100

    plot_df = pd.DataFrame(
        {
            "PC1": coords[:, pc1 - 1],
            "PC2": coords[:, pc2 - 1],
            "Sample": samples,
            "Group": [sample_groups.get(s, "Unknown") for s in samples],
        }
    )

    fig = go.Figure()
    for group, sub in plot_df.groupby("Group"):
        color = (color_mapping or {}).get(group)
        fig.add_trace(
            go.Scatter(
                x=sub["PC1"],
                y=sub["PC2"],
                mode="markers+text" if show_labels else "markers",
                name=group,
                text=sub["Sample"] if show_labels else None,
                textposition="top center",
                marker=dict(size=10, color=color),
            )
        )
        if show_ellipses:
            ellipse = _confidence_ellipse(sub[["PC1", "PC2"]].to_numpy())
            if ellipse is not None:
                fig.add_trace(
                    go.Scatter(
                        x=ellipse[0],
                        y=ellipse[1],
                        mode="lines",
                        line=dict(color=color, dash="dot"),
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )

    fig.update_layout(
        title=f"PCA Plot: PC{pc1} vs PC{pc2}",
        xaxis_title=f"PC{pc1} ({var_exp[pc1 - 1]:.1f}%)",
        yaxis_title=f"PC{pc2} ({var_exp[pc2 - 1]:.1f}%)",
        template="plotly_white",
        margin=dict(t=50, b=40),
    )
    return fig
