"""
Generates a synthetic dataset matching the column-naming conventions
the R app's build_main_data_index() expects, with a KNOWN injected
effect so you can validate the whole pipeline (filter -> impute ->
DEA) against a ground truth you control.

This unblocks development before you have a real sample file -- swap
this out for a real export once you have one, but keep this generator
around as a fast regression fixture regardless.

Usage:
    python scripts/generate_sample_dataset.py
Produces: data/sample/synthetic_proteomics.xlsx
"""
import numpy as np
import pandas as pd

N_PROTEINS = 500
N_SIGNIFICANT = 40  # first N_SIGNIFICANT proteins get a real injected fold change
GROUPS = {
    "Test": ["Test_1", "Test_2", "Test_3"],
    "Control": ["Control_1", "Control_2", "Control_3"],
}
OUT_PATH = "data/sample/synthetic_proteomics.xlsx"


def build_dataset(seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_samples = sum(len(cols) for cols in GROUPS.values())

    # Raw (non-log) abundances, lognormal-ish like real MS intensity data.
    log2_base = rng.normal(loc=20, scale=2, size=(N_PROTEINS, n_samples))
    sample_cols = [c for cols in GROUPS.values() for c in cols]
    abundance = pd.DataFrame(
        2 ** log2_base, columns=[f"Abundance: {c}" for c in sample_cols]
    )

    # Inject a real 2x (log2FC=1) increase in Test for the first N_SIGNIFICANT rows.
    test_abundance_cols = [f"Abundance: {c}" for c in GROUPS["Test"]]
    abundance.loc[: N_SIGNIFICANT - 1, test_abundance_cols] *= 2

    # Sprinkle some missingness (proteomics data is never fully observed).
    missing_mask = rng.random(abundance.shape) < 0.08
    abundance = abundance.mask(missing_mask)

    # "Found in Sample Group" columns -- boolean-ish presence flags per group,
    # matching the R app's expected header format exactly.
    group_flag_cols = {}
    for group_name, cols in GROUPS.items():
        present = abundance[[f"Abundance: {c}" for c in cols]].notna().any(axis=1)
        group_flag_cols[f"Found in Sample Group: {group_name}"] = present.map(
            {True: "High", False: "Not Found"}
        )

    metadata = pd.DataFrame(
        {
            "# Protein": range(1, N_PROTEINS + 1),
            "Accession": [f"P{10000 + i}" for i in range(N_PROTEINS)],
            "Description": [
                f"Synthetic protein {i} OS=Homo sapiens GN=GENE{i}" for i in range(N_PROTEINS)
            ],
            "Gene Symbol": [f"GENE{i}" for i in range(N_PROTEINS)],
        }
    )

    df = pd.concat(
        [metadata, pd.DataFrame(group_flag_cols), abundance],
        axis=1,
    )
    # Attach ground truth as a sidecar column for test assertions --
    # remove or move to a separate sheet if this needs to look like a
    # truly "real" export with no extra columns.
    df["_ground_truth_significant"] = df.index < N_SIGNIFICANT
    return df


if __name__ == "__main__":
    import os

    os.makedirs("data/sample", exist_ok=True)
    dataset = build_dataset()
    dataset.to_excel(OUT_PATH, index=False)
    print(f"Wrote {len(dataset)} rows to {OUT_PATH}")
    print(f"Columns: {list(dataset.columns)[:10]} ... ({len(dataset.columns)} total)")
