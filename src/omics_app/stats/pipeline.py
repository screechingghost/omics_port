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

import re

import numpy as np
import pandas as pd

from omics_app.data.columns import build_main_data_index, extract_group_names_from_columns
from omics_app.data.filtering import filter_valids, impute_downshift
from omics_app.stats.anova_path import run_anova_per_row
from omics_app.stats.mbqn import mbqn_median
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
    log2_data = pd.DataFrame(np.log2(raw), index=raw.index, columns=raw.columns)

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
    normalized = mbqn_median(imputed)

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
        # Column-name lists for the normalized abundance data, added for
        # ui/visualization.py: it needs to reconstruct R's
        # intensity_matrix (all kept proteins x samples) directly from
        # this dict's "results_df", and these are the only reliable way
        # to find those columns again (they're not distinguishable from
        # metadata/raw-abundance columns by name pattern alone).
        "test_col_names": test_col_names,
        "control_col_names": control_col_names,
    }


def _match_abundance_by_group(
    group_names: list[str], abundance_cols: list[str]
) -> dict[str, list[str]]:
    """
    Matches abundance columns to each group name via case-insensitive
    substring search, same approach R uses everywhere it auto-detects
    abundance columns from a group name (e.g. R lines 3181, 3422, 4025,
    4937 -- all plain `grep(group_name, abundance_cols, ...)`).

    That plain substring approach is ambiguous whenever one group name is
    a substring of another -- e.g. "IRRADIATED" vs "NON_IRRADIATED", or
    "CONTROL" vs "CONTROL_TREATED". A naive per-group match (each group
    searched independently against the full column list) lets the
    shorter name's pattern also match the longer name's columns, so the
    same raw column ends up claimed by two groups; when those columns are
    later renamed per-group, the second group's rename silently clobbers
    the first group's mapping for that column, corrupting column counts
    with no error until a downstream KeyError -- exactly what this
    function exists to prevent.

    Fix: process group names longest-first and remove already-claimed
    columns from the pool before matching shorter names, so a more
    specific name always wins over one that's merely a substring of it.
    This makes matching deterministic regardless of group order in the
    UI, at the cost of no longer exactly replicating R's per-group
    behavior in ambiguous cases -- R has the same underlying ambiguity,
    just resolved arbitrarily by dict/list evaluation order instead of
    by name specificity.
    """
    remaining = list(abundance_cols)
    matched_by_group: dict[str, list[str]] = {}
    for name in sorted(group_names, key=len, reverse=True):
        pattern = re.escape(name)
        matched = [c for c in remaining if re.search(pattern, c, re.IGNORECASE)]
        matched_by_group[name] = matched
        claimed = set(matched)
        remaining = [c for c in remaining if c not in claimed]
    # Return in the caller's original group order, not longest-first.
    return {name: matched_by_group[name] for name in group_names}


def run_multi_group_analysis(
    df: pd.DataFrame,
    comparison: dict,
    min_valid_percent: float,
    pvalue_threshold: float,
    significance_method: str,  # "fdr" or "raw"
    expected_replicate_counts: dict[str, int] | None = None,
) -> dict:
    """
    Port of run_anova_comparison (R lines 4012-4160) -- N-group one-way
    ANOVA per row (log2 -> filter -> impute -> normalize -> aov(), same
    shape as run_two_group_analysis above, using run_anova_per_row for
    step 5 instead of limma).

    comparison: one entry from store-comparisons' "comparisons" list,
        method == "anova". Must have "groups": the raw "Found in Sample
        Group: X" column headers selected on the Comparisons tab.

    Unlike R (which accepts a manually-specified abundance_cols_by_group
    from a per-group column picker -- R lines 4023-4033), no such picker
    is wired in ui/comparisons.py yet, so this always auto-matches
    abundance columns to each group name -- see _match_abundance_by_group
    for how that matching handles one group name being a substring of
    another (e.g. "IRRADIATED" vs "NON_IRRADIATED"), which R's own
    per-group grep() does not handle safely.

    No fold-change concept for ANOVA -- an F-test across 3+ groups has no
    single "up/down" direction, so only the p-value/q-value threshold
    applies (matches the report generator's "Not applicable for global
    ANOVA" text for this case).

    expected_replicate_counts: optional {group_name: count}, ported from
    R's "Replicate Counts" textInput + "Auto-Detect Replicates" button
    (R lines 2892-2901, 3222-3261; validated at R lines 4035-4041). Purely
    a safety check -- if given, raises AnalysisError when what's
    auto-matched right now doesn't match what the user expects, instead
    of silently running on however many columns happened to match (catches
    a typo'd group name, or a sample that dropped out of the export
    leaving fewer/extra columns than intended). Omit it to keep the old
    behavior: whatever auto-matches is used with no check.

    Returns a dict: results_df, n_kept, n_significant(_pvalue/_qvalue),
    group_names, sig_column_used, group_col_names (renamed normalized-
    abundance columns per group, the ANOVA analogue of test_col_names/
    control_col_names -- see run_two_group_analysis for why renaming is
    needed at all).
    """
    group_cols = comparison.get("groups") or []
    if len(group_cols) < 2:
        raise AnalysisError("ANOVA requires at least 2 groups.")

    group_names = extract_group_names_from_columns(group_cols)
    if len(group_names) < 2:
        raise AnalysisError("ANOVA requires at least 2 distinct groups.")

    abundance_cols_all = build_main_data_index(df)["abundance_cols"]
    abundance_by_group = _match_abundance_by_group(group_names, abundance_cols_all)
    empty_groups = [g for g, cols in abundance_by_group.items() if not cols]
    if empty_groups:
        raise AnalysisError(
            f"No abundance columns auto-detected for group(s): {', '.join(empty_groups)}."
        )

    if expected_replicate_counts:
        mismatches = [
            f"{g}: expected {expected_replicate_counts[g]}, found {len(abundance_by_group[g])}"
            for g in group_names
            if g in expected_replicate_counts
            and expected_replicate_counts[g] != len(abundance_by_group[g])
        ]
        if mismatches:
            raise AnalysisError(
                "Mismatch between replicate counts and detected columns -- "
                + "; ".join(mismatches)
                + ". Check your group names, or click Auto-Detect Replicates again."
            )

    all_abundance = [c for cols in abundance_by_group.values() for c in cols]
    replicate_counts = {g: len(cols) for g, cols in abundance_by_group.items()}

    missing = [c for c in all_abundance if c not in df.columns]
    if missing:
        raise AnalysisError(f"Columns not found in dataset: {missing}")

    # 1. log2 transform
    raw = df[all_abundance].copy()
    raw = raw.mask(raw <= 0)
    log2_data = np.log2(raw)

    # 2. Filter
    min_ratio = min_valid_percent / 100
    min_count = {g: max(1, int(np.floor(n * min_ratio))) for g, n in replicate_counts.items()}
    keep_mask = filter_valids(log2_data, abundance_by_group, min_count, at_least_one=False)
    log2_filtered = log2_data.loc[keep_mask]
    n_kept = len(log2_filtered)
    if n_kept == 0:
        raise AnalysisError("No proteins passed filtering!")

    # 3. Imputation -- same fixed seed rationale as run_two_group_analysis.
    imputed = impute_downshift(log2_filtered, abundance_by_group, random_state=1)

    # 4. Normalization -- see module docstring re: mbqn() gap
    normalized = mbqn_median(imputed)

    # Rename to distinct, group-labeled column names before merging --
    # same duplicate-column-name hazard as the 2-group path.
    group_col_names: dict[str, list[str]] = {}
    rename_map: dict[str, str] = {}
    for g, cols in abundance_by_group.items():
        renamed = [f"{g}_{i + 1}" for i in range(len(cols))]
        group_col_names[g] = renamed
        rename_map.update(dict(zip(cols, renamed)))
    normalized = normalized.rename(columns=rename_map)

    # 5. ANOVA (R lines 4084-4096: aov(protein_data ~ group_factor) per row)
    ordered_cols = [c for g in group_names for c in group_col_names[g]]
    group_labels = [g for g in group_names for _ in range(replicate_counts[g])]
    anova_results = run_anova_per_row(normalized[ordered_cols], group_labels)

    # 6. Merge with metadata
    metadata_cols = [c for c in df.columns if c not in all_abundance]
    merged = pd.concat(
        [
            df.loc[anova_results.index, metadata_cols],
            anova_results[["F_statistic", "P.Value", "adj.P.Val"]],
            normalized.loc[anova_results.index],
            df.loc[anova_results.index, all_abundance],
        ],
        axis=1,
    )

    # 7. Significance -- no fold-change filter option for ANOVA (R line
    # 4117-4119 doesn't apply one either).
    merged["Significance_pvalue"] = merged["P.Value"] < pvalue_threshold
    merged["Significance_qvalue"] = merged["adj.P.Val"] < pvalue_threshold
    merged["Groups"] = "_vs_".join(group_names)
    merged["Comparison"] = f"ANOVA: {' vs '.join(group_names)}"

    sig_column = "Significance_qvalue" if significance_method == "fdr" else "Significance_pvalue"
    n_significant = int(merged[sig_column].sum())
    n_significant_pvalue = int(merged["Significance_pvalue"].sum())
    n_significant_qvalue = int(merged["Significance_qvalue"].sum())

    return {
        "results_df": merged,
        "n_kept": n_kept,
        "n_significant": n_significant,
        "n_significant_pvalue": n_significant_pvalue,
        "n_significant_qvalue": n_significant_qvalue,
        "group_names": group_names,
        "sig_column": sig_column,
        # ANOVA analogue of test_col_names/control_col_names -- one
        # renamed-column list per group instead of exactly two.
        "group_col_names": group_col_names,
    }
