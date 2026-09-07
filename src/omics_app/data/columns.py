"""
Port of app_12-02.R functions:
  extract_group_names_from_columns  (line 69)
  build_main_data_index             (line 1798)

This is the auto-detection logic that lets the app work generically
across different Proteome-Discoverer-style exports: it classifies each
column header into group / abundance / metadata based on naming
patterns, without the user having to hand-map every column.
"""
import re

import pandas as pd


def build_main_data_index(df: pd.DataFrame) -> dict:
    """
    Classifies columns of an uploaded data file into:
      group_cols           -- "Found in Sample Group: X" columns
      abundance_cols       -- "Abundance...: X" columns (excluding
                               Ratio/Normalized/Scaled/Grouped variants)
      metadata_candidates  -- everything else (Accession, Description,
                               Gene Symbol, etc.)
      default_selected     -- same as metadata_candidates in the R app

    Mirrors build_main_data_index exactly, including its exclude-pattern
    order (this matters: "Abundance" is excluded from metadata_candidates
    even for columns that didn't match the abundance_cols regex, e.g. a
    literal column named "Abundance Count").
    """
    if df is None:
        return {
            "group_cols": [],
            "abundance_cols": [],
            "metadata_candidates": [],
            "default_selected": [],
        }

    all_cols = list(df.columns)

    group_cols = [c for c in all_cols if re.search(r"Found in Sample.*Group", c, re.IGNORECASE)]

    abundance_cols = [c for c in all_cols if re.search(r"Abundance.*:", c, re.IGNORECASE)]
    abundance_cols = [
        c for c in abundance_cols
        if not re.search(r"Ratio|Normalized|Scaled|Grouped", c, re.IGNORECASE)
    ]

    exclude_patterns = ["Abundance", "Found in Sample", "# Protein"]
    metadata_candidates = list(all_cols)
    for pattern in exclude_patterns:
        metadata_candidates = [
            c for c in metadata_candidates if not re.search(pattern, c, re.IGNORECASE)
        ]

    return {
        "group_cols": group_cols,
        "abundance_cols": abundance_cols,
        "metadata_candidates": metadata_candidates,
        "default_selected": list(metadata_candidates),
    }


def extract_group_names_from_columns(group_columns: list[str]) -> list[str]:
    """
    Port of extract_group_names_from_columns (line 69). Pulls the group
    name out of "Found in Sample Group: GROUPNAME" style headers; falls
    back to the raw column name if it doesn't match the pattern.
    """
    names: list[str] = []
    for col in group_columns:
        if re.search(r"Found in Sample.*Group.*:", col, re.IGNORECASE):
            parts = col.split(":", 1)
            if len(parts) == 2:
                names.append(parts[1].strip())
                continue
        names.append(col)
    # dedupe while preserving order (R's unique() keeps first-seen order)
    seen = set()
    unique_names = []
    for n in names:
        if n not in seen:
            seen.add(n)
            unique_names.append(n)
    return unique_names
