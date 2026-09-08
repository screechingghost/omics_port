import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def run_anova_per_row(
    intensity_matrix: pd.DataFrame,
    group_labels: list[str],
) -> pd.DataFrame:
    """
    intensity_matrix: rows = proteins, columns = samples (already log2,
                       filtered, imputed, normalized -- upstream of this).
    group_labels: group name for each column, same order/length as columns.

    Returns a DataFrame indexed like intensity_matrix with F_statistic,
    P.Value, adj.P.Val (BH), and Mean_<group> for each group.
    """
    groups = pd.Series(group_labels)
    unique_groups = groups.unique()

    f_stats = np.empty(len(intensity_matrix))
    p_values = np.empty(len(intensity_matrix))

    values = intensity_matrix.to_numpy()
    for i in range(values.shape[0]):
        row = values[i]
        samples_by_group = [row[groups.values == g] for g in unique_groups]
        f_stat, p_val = stats.f_oneway(*samples_by_group)
        f_stats[i] = f_stat
        p_values[i] = p_val

    adj_p = multipletests(p_values, method="fdr_bh")[1]

    result = pd.DataFrame(
        {"F_statistic": f_stats, "P.Value": p_values, "adj.P.Val": adj_p},
        index=intensity_matrix.index,
    )
    for g in unique_groups:
        cols = intensity_matrix.columns[groups.values == g]
        result[f"Mean_{g}"] = intensity_matrix[cols].mean(axis=1)
    return result
