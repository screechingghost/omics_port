import numpy as np
import pandas as pd

from omics_app.data.filtering import filter_valids, parse_column_input


def test_parse_column_input_ranges_and_singles():
    assert parse_column_input("1-5,8,10") == [1, 2, 3, 4, 5, 8, 10]


def test_parse_column_input_single_value():
    assert parse_column_input("7") == [7]


def test_filter_valids_requires_all_groups_by_default():
    df = pd.DataFrame(
        {
            "test_1": [1.0, np.nan],
            "test_2": [1.0, np.nan],
            "control_1": [1.0, 1.0],
            "control_2": [1.0, 1.0],
        }
    )
    group_columns = {"test": ["test_1", "test_2"], "control": ["control_1", "control_2"]}
    keep = filter_valids(df, group_columns, min_count={"test": 2, "control": 2})
    assert keep.tolist() == [True, False]


def test_filter_valids_at_least_one_group():
    df = pd.DataFrame(
        {
            "test_1": [np.nan],
            "test_2": [np.nan],
            "control_1": [1.0],
            "control_2": [1.0],
        }
    )
    group_columns = {"test": ["test_1", "test_2"], "control": ["control_1", "control_2"]}
    keep = filter_valids(
        df, group_columns, min_count={"test": 2, "control": 2}, at_least_one=True
    )
    assert keep.tolist() == [True]
