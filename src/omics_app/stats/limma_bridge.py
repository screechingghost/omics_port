"""
Port of run_normal_comparison (app_12-02.R, line 3845) -- the 2-group
path. This is the ONLY part of the stats engine that needs limma; keep
it running in R via rpy2 rather than reimplementing eBayes moderation.

Requires: R installed with the `limma` Bioconductor package, and rpy2
installed in this venv (see README.md "R + rpy2 setup").
"""
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()
limma = importr("limma")
base = importr("base")
stats_r = importr("stats")


def run_limma_two_group(
    intensity_matrix: pd.DataFrame,
    test_group: str,
    control_group: str,
    group_labels: list[str],
) -> pd.DataFrame:
    """
    intensity_matrix: rows = proteins, columns = samples (log2, filtered,
                       imputed, normalized -- matches R's pd_output_fin).
    group_labels: "test" / "control" for each column, same order as columns.

    Returns a DataFrame indexed like intensity_matrix with logFC,
    P.Value, adj.P.Val -- equivalent to R's topTable() output.
    """
    with (robjects.default_converter + pandas2ri.converter).context():
        r_matrix = robjects.conversion.py2rpy(intensity_matrix)

    factor_r = robjects.FactorVector(
        robjects.StrVector(group_labels), levels=robjects.StrVector(["test", "control"])
    )
    design = stats_r.model_matrix(robjects.Formula("~0+factor_r"), **{"factor_r": factor_r})
    design.colnames = robjects.StrVector([test_group, control_group])

    contrast_str = f"`{test_group}`-`{control_group}`"
    contrast_matrix = limma.makeContrasts(contrasts=contrast_str, levels=design)

    fit = limma.lmFit(r_matrix, design)
    fit_contrasts = limma.contrasts_fit(fit, contrast_matrix)
    fit_eb = limma.eBayes(fit_contrasts)

    top_table = limma.topTable(fit_eb, number=intensity_matrix.shape[0])
    with (robjects.default_converter + pandas2ri.converter).context():
        result = robjects.conversion.rpy2py(top_table)

    return result[["logFC", "P.Value", "adj.P.Val"]]
