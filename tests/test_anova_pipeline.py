import numpy as np
import pandas as pd
import pytest

from omics_app.stats.pipeline import AnalysisError, run_multi_group_analysis


def _synthetic_dataset(seed: int = 0, n_proteins: int = 60, n_per_group: int = 4):
    rng = np.random.RandomState(seed)
    groups = ["GroupA", "GroupB", "GroupC"]
    data = {}
    for g in groups:
        for i in range(1, n_per_group + 1):
            data[f"Abundance: F{i}, Sample, {g}"] = rng.lognormal(
                mean=10, sigma=1, size=n_proteins
            )
    df = pd.DataFrame(data)
    df["Accession"] = [f"P{i:04d}" for i in range(n_proteins)]
    df["Description"] = [f"protein {i}" for i in range(n_proteins)]
    df = df.set_index("Accession")
    return df, groups


def _comparison(groups):
    return {
        "name": "Comp1",
        "method": "anova",
        "groups": [f"Found in Sample Group: {g}" for g in groups],
    }


def test_requires_at_least_two_groups():
    df, _groups = _synthetic_dataset()
    comparison = {"groups": ["Found in Sample Group: GroupA"]}
    with pytest.raises(AnalysisError, match="at least 2 groups"):
        run_multi_group_analysis(df, comparison, 70, 0.05, "raw")


def test_missing_abundance_columns_for_a_group_raises():
    df, _groups = _synthetic_dataset()
    comparison = {
        "groups": [
            "Found in Sample Group: GroupA",
            "Found in Sample Group: DoesNotExist",
        ]
    }
    with pytest.raises(AnalysisError, match="No abundance columns auto-detected"):
        run_multi_group_analysis(df, comparison, 70, 0.05, "raw")


def test_runs_and_produces_expected_shape():
    df, groups = _synthetic_dataset()
    result = run_multi_group_analysis(
        df, _comparison(groups), min_valid_percent=70, pvalue_threshold=0.05,
        significance_method="raw",
    )
    assert result["group_names"] == groups
    assert set(result["group_col_names"]) == set(groups)
    for g in groups:
        assert result["group_col_names"][g] == [f"{g}_1", f"{g}_2", f"{g}_3", f"{g}_4"]

    results_df = result["results_df"]
    for col in ["F_statistic", "P.Value", "adj.P.Val", "Significance_pvalue",
                "Significance_qvalue", "Groups", "Comparison"]:
        assert col in results_df.columns
    assert (results_df["Groups"] == "GroupA_vs_GroupB_vs_GroupC").all()
    assert (results_df["Comparison"] == "ANOVA: GroupA vs GroupB vs GroupC").all()
    assert result["n_kept"] == len(results_df)
    assert result["sig_column"] == "Significance_pvalue"


def test_significance_method_selects_sig_column():
    df, groups = _synthetic_dataset()
    result = run_multi_group_analysis(
        df, _comparison(groups), min_valid_percent=70, pvalue_threshold=0.05,
        significance_method="fdr",
    )
    assert result["sig_column"] == "Significance_qvalue"


def test_detects_injected_group_difference():
    df, groups = _synthetic_dataset(n_proteins=100)
    signal_rows = df.index[:15]
    group_a_cols = [c for c in df.columns if "GroupA" in c]
    df.loc[signal_rows, group_a_cols] = df.loc[signal_rows, group_a_cols] * 6

    result = run_multi_group_analysis(
        df, _comparison(groups), min_valid_percent=70, pvalue_threshold=0.05,
        significance_method="raw",
    )
    sig_df = result["results_df"]
    caught = sig_df.index[sig_df["Significance_pvalue"]]
    # Not every injected row needs to clear the threshold with only 4
    # replicates/group, but the large injected effect should recover a
    # clear majority of them -- a regression that broke row alignment
    # between filter/impute/normalize/ANOVA steps would make this ~0.
    overlap = len(set(signal_rows) & set(caught))
    assert overlap >= len(signal_rows) * 0.5


def test_group_name_substring_of_another_does_not_cross_contaminate():
    """Regression test: a group name that's a substring of another group's
    name (e.g. "IRRADIATED" vs "NON_IRRADIATED") must not have its
    abundance columns stolen by/merged with the other group. Previously,
    naive per-group substring matching let "IRRADIATED"'s pattern also
    match "NON_IRRADIATED" columns, corrupting column counts and raising
    a KeyError once the ANOVA step tried to select the (silently
    clobbered) renamed columns."""
    rng = np.random.RandomState(1)
    n_proteins = 50
    groups = {"IRRADIATED": 3, "NON_IRRADIATED": 3}
    data = {}
    for g, n in groups.items():
        for i in range(1, n + 1):
            data[f"Abundance: F{i}, Sample, {g}"] = rng.lognormal(10, 1, size=n_proteins)
    df = pd.DataFrame(data)
    df["Accession"] = [f"P{i:04d}" for i in range(n_proteins)]
    df = df.set_index("Accession")

    comparison = {"groups": [f"Found in Sample Group: {g}" for g in groups]}
    result = run_multi_group_analysis(
        df, comparison, min_valid_percent=70, pvalue_threshold=0.05, significance_method="raw"
    )

    assert result["group_col_names"]["IRRADIATED"] == [
        "IRRADIATED_1", "IRRADIATED_2", "IRRADIATED_3"
    ]
    assert result["group_col_names"]["NON_IRRADIATED"] == [
        "NON_IRRADIATED_1", "NON_IRRADIATED_2", "NON_IRRADIATED_3"
    ]
    # Every ordered column referenced by the ANOVA step must actually
    # exist in the results -- this is what raised KeyError before the fix.
    ordered_cols = (
        result["group_col_names"]["IRRADIATED"] + result["group_col_names"]["NON_IRRADIATED"]
    )
    for col in ordered_cols:
        assert col in result["results_df"].columns


def test_expected_replicate_counts_pass_when_matching():
    df, groups = _synthetic_dataset()
    result = run_multi_group_analysis(
        df, _comparison(groups), min_valid_percent=70, pvalue_threshold=0.05,
        significance_method="raw",
        expected_replicate_counts={g: 4 for g in groups},
    )
    assert result["n_kept"] > 0


def test_expected_replicate_counts_raise_on_mismatch():
    df, groups = _synthetic_dataset()
    with pytest.raises(AnalysisError, match="Mismatch between replicate counts"):
        run_multi_group_analysis(
            df, _comparison(groups), min_valid_percent=70, pvalue_threshold=0.05,
            significance_method="raw",
            expected_replicate_counts={"GroupA": 5, "GroupB": 4, "GroupC": 4},
        )
