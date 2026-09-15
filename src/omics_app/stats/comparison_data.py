"""
Shared helper for reconstructing one comparison's data from the flat
dict ui/analysis.py's run_analysis callback stores in
store-dea-results. Used by ui/enrichment.py; ui/visualization.py has
its own local copy of this same logic (_load_comparison) rather than
importing this module -- consolidating that is a small future
cleanup, not required for either tab to work correctly.
"""

import pandas as pd


def load_comparison(dea_results: dict | None, comp_name: str | None):
    """
    Reconstructs (df, sample_columns, sample_groups, sig_mask,
    test_group_name, control_group_name) for one comparison. Returns
    None if the comparison is missing, or isn't the 2-group ("normal")
    method this port supports.
    """
    if not dea_results or not comp_name or comp_name not in dea_results:
        return None
    entry = dea_results[comp_name]
    if "test_col_names" not in entry or "control_col_names" not in entry:
        # Either an ANOVA/uploaded-only entry, or from before
        # pipeline.py/analysis.py were extended to store these --
        # re-run Analysis for this comparison to populate them.
        return None

    df = pd.DataFrame(entry["results"])
    rowname_col = entry.get("rowname_col")
    if rowname_col and rowname_col in df.columns:
        df = df.set_index(rowname_col)

    test_cols = entry["test_col_names"]
    control_cols = entry["control_col_names"]
    sample_groups = {c: entry["test_group_name"] for c in test_cols}
    sample_groups.update({c: entry["control_group_name"] for c in control_cols})

    sig_col = entry.get("sig_column", "Significance_pvalue")
    sig_mask = df[sig_col] if sig_col in df.columns else pd.Series(False, index=df.index)

    return {
        "df": df,
        "sample_columns": test_cols + control_cols,
        "sample_groups": sample_groups,
        "sig_mask": sig_mask,
        "test_group_name": entry["test_group_name"],
        "control_group_name": entry["control_group_name"],
    }
