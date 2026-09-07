"""
Port of observeEvent(input$run_analysis_btn, ...) (app_12-02.R, lines
4339-4600) -- the 2-group ("normal") analysis pipeline: log2 transform
-> filter -> impute -> normalize -> limma DEA -> merge -> significance.

Deliberately kept UI-free and pure (DataFrame + dict in, dict out) so it
can be unit/integration tested directly, same as every other stats
module in this project -- ui/analysis.py is a thin wrapper around this.

Known gap: R normalizes with mbqn() (Median Balanced Quantile
Normalization, from the Bioconductor/CRAN `MBQN` package). No direct
Python equivalent was found, so this uses simple per-sample median
centering as a placeholder (shift each sample so its median matches the
grand median). This is NOT numerically equivalent to MBQN -- if you
need exact parity with the R app's output, porting the actual MBQN
algorithm is the next thing to check, not something to assume this
placeholder covers.
"""

import numpy as np
import pandas as pd

from omics_app.data.columns import extract_group_names_from_columns
from omics_app.data.filtering import filter_valids, impute_downshift
from omics_app.stats.pylimma_bridge import run_pylimma_two_group


class AnalysisError(Exception):
    """Raised for user-facing validation failures (missing columns,
    nothing passed filtering, etc.) -- caught by the UI layer and shown
    as an alert, mirroring R's showNotification(..., type="error")."""


def _median_center_normalize(log2_data: pd.DataFrame) -> pd.DataFrame:
    """Placeholder for R's mbqn() -- see module docstring."""
    grand_median = log2_data.stack().median()
    return log2_data.apply(lambda col: col - col.median() + grand_median, axis=0)


def run_two_group_analysis(
    df: pd.DataFrame,
    comparison: dict,
    min_valid_percent: float,
    pvalue_threshold: float,
    significance_method: str,  # "fdr" or "raw"
    fc_threshold: float = 0,
    fc_operator: str = "none",  # none/gte/gt/lte/lt/abs
) -> dict:
    """
    df: the full main dataset, indexed by protein ID (rowname_col already
        set as index -- caller's responsibility, matching how
        store-main-data's rowname_col is used elsewhere).
    comparison: one entry from store-comparisons' "comparisons" list,
        method == "normal" (2-group). Must have group1/group2 (the raw
        "Found in Sample Group: X" column names) and abundance_g1/
        abundance_g2 (the actual abundance columns for each group) --
        both already collected during Comparisons tab Setup, so this
        does NOT re-run the auto-detection R does inside the Analysis
        tab itself (see ui/analysis.py module docstring).

    Returns a dict: results_df, n_kept, n_significant, test_group_name,
    control_group_name, sig_column_used.
    """
    abundance_test = comparison.get("abundance_g1") or []
    abundance_control = comparison.get("abundance_g2") or []
    if not abundance_test:
        raise AnalysisError("No test-group abundance columns selected for this comparison.")
    if not abundance_control:
        raise AnalysisError("No control-group abundance columns selected for this comparison.")

    test_group_name = extract_group_names_from_columns([comparison["group1"]])[0]
    control_group_name = extract_group_names_from_columns([comparison["group2"]])[0]

    n_test, n_control = len(abundance_test), len(abundance_control)
    all_abundance = abundance_test + abundance_control

    missing = [c for c in all_abundance if c not in df.columns]
    if missing:
        raise AnalysisError(f"Columns not found in dataset: {missing}")

    # 1. log2 transform (R lines 4404-4410): non-positive values -> NaN first.
    raw = df[all_abundance].copy()
    raw = raw.mask(raw <= 0)
    log2_data = np.log2(raw)

    # 2. Filter (R lines 4412-4436)
    min_ratio = min_valid_percent / 100
    group_columns = {test_group_name: abundance_test, control_group_name: abundance_control}
    min_count = {
        test_group_name: max(1, int(np.floor(n_test * min_ratio))),
        control_group_name: max(1, int(np.floor(n_control * min_ratio))),
    }
    keep_mask = filter_valids(log2_data, group_columns, min_count, at_least_one=False)
    log2_filtered = log2_data.loc[keep_mask]
    n_kept = len(log2_filtered)
    if n_kept == 0:
        raise AnalysisError("No proteins passed filtering!")

    # 3. Imputation (R lines 4453-4458)
    # random_state=1 matches R's set.seed(1) (called before every
    # imputation, R lines 135/156) -- without a fixed seed here, every
    # click of "Run Analysis" would silently produce different imputed
    # values and therefore different significant-protein counts even
    # with identical settings, which is not what R does and would look
    # like a bug (results changing for no visible reason).
    imputed = impute_downshift(log2_filtered, group_columns, random_state=1)

    # 4. Normalization -- see module docstring re: mbqn() gap
    normalized = _median_center_normalize(imputed)

    # Rename normalized columns to distinct names before merging, matching
    # R's dummy_names approach (lines 4395-4402): without this, the
    # normalized columns and the raw abundance columns appended later
    # would have IDENTICAL names (both start from the same abundance
    # column names), silently corrupting to_dict()/to_records() calls
    # downstream (pandas drops data on duplicate column names).
    test_col_names = [f"{test_group_name}_{i + 1}" for i in range(n_test)]
    control_col_names = [f"{control_group_name}_{i + 1}" for i in range(n_control)]
    normalized = normalized.rename(
        columns=dict(zip(abundance_test + abundance_control, test_col_names + control_col_names))
    )

    # 5. DEA (R lines 4466-4488)
    group_labels = ["test"] * n_test + ["control"] * n_control
    dea = run_pylimma_two_group(
        intensity_matrix=normalized[test_col_names + control_col_names],
        test_group=test_group_name,
        control_group=control_group_name,
        group_labels=group_labels,
    )

    # 6. Merge with metadata (R lines 4492-4506)
    metadata_cols = [c for c in df.columns if c not in all_abundance]
    merged = pd.concat(
        [
            df.loc[dea.index, metadata_cols],
            dea[["logFC", "P.Value", "adj.P.Val"]],
            normalized.loc[dea.index],
            df.loc[dea.index, all_abundance],
        ],
        axis=1,
    )

    # 7. Significance (R lines 4508-4528)
    pval_sig = merged["P.Value"] < pvalue_threshold
    qval_sig = merged["adj.P.Val"] < pvalue_threshold

    if fc_operator != "none" and fc_threshold > 0:
        fc_ops = {
            "gte": lambda x: x >= fc_threshold,
            "gt": lambda x: x > fc_threshold,
            "lte": lambda x: x <= -fc_threshold,
            "lt": lambda x: x < -fc_threshold,
            "abs": lambda x: x.abs() >= fc_threshold,
        }
        fc_sig = fc_ops.get(fc_operator, lambda x: pd.Series(True, index=x.index))(merged["logFC"])
        merged["Significance_pvalue"] = pval_sig & fc_sig
        merged["Significance_qvalue"] = qval_sig & fc_sig
    else:
        merged["Significance_pvalue"] = pval_sig
        merged["Significance_qvalue"] = qval_sig

    merged["Test_Group"] = test_group_name
    merged["Control_Group"] = control_group_name
    merged["Comparison"] = f"{test_group_name} vs {control_group_name}"

    sig_column = "Significance_qvalue" if significance_method == "fdr" else "Significance_pvalue"
    n_significant = int(merged[sig_column].sum())

    # R's summary (lines 5370-5378) always shows BOTH counts regardless
    # of which one is used for the table filter, plus up/down regulation
    # counts computed specifically from the p-value-significant set
    # (not whichever sig_column was chosen) -- matched exactly here.
    n_significant_pvalue = int(merged["Significance_pvalue"].sum())
    n_significant_qvalue = int(merged["Significance_qvalue"].sum())
    pvalue_sig_rows = merged[merged["Significance_pvalue"]]
    n_upregulated = int((pvalue_sig_rows["logFC"] > 0).sum())
    n_downregulated = int((pvalue_sig_rows["logFC"] < 0).sum())

    return {
        "results_df": merged,
        "n_kept": n_kept,
        "n_significant": n_significant,
        "n_significant_pvalue": n_significant_pvalue,
        "n_significant_qvalue": n_significant_qvalue,
        "n_upregulated": n_upregulated,
        "n_downregulated": n_downregulated,
        "test_group_name": test_group_name,
        "control_group_name": control_group_name,
        "sig_column": sig_column,
    }
