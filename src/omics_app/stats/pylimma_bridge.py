"""
Port of run_normal_comparison limma calls (app_12-02.R, lines 3893-3903)
using pylimma: no-intercept two-group design, contrast test - control,
lmFit -> contrasts.fit -> eBayes -> topTable.
"""
import numpy as np
import pandas as pd
from pylimma import contrasts_fit, e_bayes, lm_fit, make_contrasts, top_table

_COLUMN_RENAME = {
    "log_fc": "logFC",
    "p_value": "P.Value",
    "adj_p_value": "adj.P.Val",
}


def run_pylimma_two_group(
    intensity_matrix: pd.DataFrame,
    test_group: str,
    control_group: str,
    group_labels: list[str],
) -> pd.DataFrame:
    """
    intensity_matrix: rows = proteins, columns = samples (log2, filtered,
                       imputed, normalized -- matches R's pd_output_fin).
    group_labels: "test" / "control" for each column, same order as
                  intensity_matrix.columns.

    Returns a DataFrame indexed like intensity_matrix with logFC,
    P.Value, adj.P.Val -- renamed to match R topTable() conventions.
    """
    labels = np.asarray(group_labels)
    design = np.column_stack(
        [
            (labels == "test").astype(float),
            (labels == "control").astype(float),
        ]
    )
    contrast = make_contrasts(
        f"{test_group}-{control_group}",
        levels=[test_group, control_group],
    )

    fit = lm_fit(intensity_matrix.to_numpy(dtype=float), design)
    fit = contrasts_fit(fit, contrasts=contrast.to_numpy())
    fit = e_bayes(fit)

    result = top_table(
        fit,
        coef=0,
        number=intensity_matrix.shape[0],
        sort_by="none",
    )
    result = result.rename(columns=_COLUMN_RENAME)
    result.index = intensity_matrix.index
    return result[["logFC", "P.Value", "adj.P.Val"]]
