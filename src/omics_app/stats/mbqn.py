"""
Port of MBQN::mbqn(x, FUN = median, method = "limma") used in
app_12-02.R (run_normal_comparison, line 3890).

Algorithm (R/mbqn.R): subtract a per-row offset (median), quantile-normalize
the residuals with limma::normalizeQuantiles, then add the offset back.
"""
import numpy as np
import pandas as pd
from pylimma.normalize import normalize_quantiles


def mbqn_median(data: pd.DataFrame) -> pd.DataFrame:
    """Median-balanced quantile normalization; rows = features, columns = samples."""
    matrix = data.to_numpy(dtype=float, copy=True)
    matrix[~np.isfinite(matrix)] = np.nan
    row_offset = np.nanmedian(matrix, axis=1)
    centered = matrix - row_offset[:, np.newaxis]
    normalized = normalize_quantiles(centered) + row_offset[:, np.newaxis]
    return pd.DataFrame(normalized, index=data.index, columns=data.columns)
