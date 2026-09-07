"""
Run this FIRST, before writing anything else in the stats module.
Confirms R, limma, and rpy2 are correctly wired together, and gives you
a known-answer sanity check on synthetic data.

Usage:
    python scripts/validate_limma_bridge.py
"""
import numpy as np
import pandas as pd

from omics_app.stats.limma_bridge import run_limma_two_group

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n_proteins = 200
    test_cols = ["test_1", "test_2", "test_3"]
    control_cols = ["control_1", "control_2", "control_3"]

    # First 20 rows: real signal (shifted mean). Rest: noise only.
    data = pd.DataFrame(
        rng.normal(10, 1, size=(n_proteins, 6)),
        columns=test_cols + control_cols,
        index=[f"PROT_{i}" for i in range(n_proteins)],
    )
    data.loc[data.index[:20], test_cols] += 3  # inject real fold change

    result = run_limma_two_group(
        intensity_matrix=data,
        test_group="test",
        control_group="control",
        group_labels=["test"] * 3 + ["control"] * 3,
    )

    print(result.head(10))
    n_sig = (result["P.Value"] < 0.05).sum()
    print(f"\n{n_sig} of {n_proteins} proteins significant at p<0.05 "
          f"(expect roughly 20 real + a handful of false positives)")
