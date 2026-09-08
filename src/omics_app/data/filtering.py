"""
Port of R functions (app_12-02.R):
  filter_valids_custom   (line 89)
  filter_valids_anova    (line 109)
  impute_data_custom     (line 129)
  impute_data_anova      (line 150)
  parse_column_input     (line 171)
"""

import numpy as np
import pandas as pd


def filter_valids(
    log2_data: pd.DataFrame,
    group_columns: dict[str, list[str]],
    min_count: dict[str, int],
    at_least_one: bool = False,
) -> pd.Series:
    """
    Returns a boolean Series (index-aligned to log2_data) marking which
    rows have enough non-missing values per group.

    group_columns: {group_name: [column names in log2_data for that group]}
    min_count: {group_name: minimum non-NA values required in that group}
    at_least_one: if True, keep row if ANY group meets min_count;
                  if False, ALL groups must meet min_count (matches R default).
    """
    keep_per_group = pd.DataFrame(
        {
            group: log2_data[cols].notna().sum(axis=1) >= min_count[group]
            for group, cols in group_columns.items()
        }
    )
    return keep_per_group.any(axis=1) if at_least_one else keep_per_group.all(axis=1)


def impute_downshift(
    log2_data: pd.DataFrame,
    group_columns: dict[str, list[str]] | None = None,
    downshift: float = 1.8,
    width: float = 0.3,
    random_state: int | None = None,
    keep: pd.Series | None = None,
) -> pd.DataFrame:
    """
    Port of impute_data_custom / impute_data_anova (R lines 129-168).

    Per *column* (sample), missing/non-finite values are drawn from
    N(mean - downshift*sd, width*sd), where mean/sd are computed from
    finite values in that column among KEEP rows. After filtering, KEEP
    is all True, so this is column-wide MNAR imputation.

    R's test_pattern / control_pattern / group_patterns arguments are
    unused in the R source; group_columns is kept only for call-site
    compatibility and is ignored.

    NumPy and R do not share an RNG, so seeds will not be bit-identical.
    """
    del group_columns  # unused, matching R
    rng = np.random.RandomState(random_state)
    result = log2_data.copy()
    values = result.to_numpy(dtype=float, copy=True)
    values[~np.isfinite(values)] = np.nan

    if keep is None:
        keep_mask = np.ones(values.shape[0], dtype=bool)
    else:
        keep_mask = keep.reindex(result.index).fillna(False).to_numpy(dtype=bool)

    for j in range(values.shape[1]):
        col = values[:, j]
        kept = col[keep_mask]
        finite = kept[np.isfinite(kept)]
        if finite.size == 0:
            continue
        col_sd = float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0
        col_mean = float(np.mean(finite))
        mu = col_mean - downshift * col_sd
        sigma = width * col_sd
        missing = ~np.isfinite(col)
        n_missing = int(missing.sum())
        if n_missing == 0:
            continue
        values[missing, j] = rng.normal(mu, sigma, size=n_missing)

    return pd.DataFrame(values, index=result.index, columns=result.columns)


def parse_column_input(text: str) -> list[int]:
    """
    Port of parse_column_input (line 171): parses strings like
    "1-5,8,10" into [1,2,3,4,5,8,10].
    """
    indices: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-")
            indices.extend(range(int(start), int(end) + 1))
        else:
            indices.append(int(part))
    return indices
