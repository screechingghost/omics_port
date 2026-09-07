import pandas as pd

from omics_app.data.columns import build_main_data_index, extract_group_names_from_columns


def _sample_columns() -> list[str]:
    return [
        "# Protein",
        "Accession",
        "Description",
        "Gene Symbol",
        "Found in Sample Group: Test",
        "Found in Sample Group: Control",
        "Abundance: Test_1",
        "Abundance: Test_2",
        "Abundance: Test_3",
        "Abundance: Control_1",
        "Abundance: Control_2",
        "Abundance: Control_3",
        "Abundance Ratio: (Test) / (Control)",  # must be excluded (Ratio)
        "Abundances Normalized: Test_1",  # must be excluded (Normalized)
    ]


def test_build_main_data_index_classifies_columns():
    df = pd.DataFrame(columns=_sample_columns())
    index = build_main_data_index(df)

    assert index["group_cols"] == [
        "Found in Sample Group: Test",
        "Found in Sample Group: Control",
    ]
    assert index["abundance_cols"] == [
        "Abundance: Test_1",
        "Abundance: Test_2",
        "Abundance: Test_3",
        "Abundance: Control_1",
        "Abundance: Control_2",
        "Abundance: Control_3",
    ]
    # Ratio/Normalized variants excluded from abundance_cols
    assert "Abundance Ratio: (Test) / (Control)" not in index["abundance_cols"]
    assert "Abundances Normalized: Test_1" not in index["abundance_cols"]

    # metadata_candidates excludes anything matching Abundance / Found in Sample / # Protein
    assert index["metadata_candidates"] == ["Accession", "Description", "Gene Symbol"]


def test_build_main_data_index_handles_none():
    index = build_main_data_index(None)
    assert index["group_cols"] == []
    assert index["metadata_candidates"] == []


def test_extract_group_names_from_columns():
    cols = ["Found in Sample Group: Test", "Found in Sample Group: Control"]
    assert extract_group_names_from_columns(cols) == ["Test", "Control"]


def test_extract_group_names_falls_back_to_raw_column():
    assert extract_group_names_from_columns(["Unrelated Column"]) == ["Unrelated Column"]
