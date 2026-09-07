"""
Run this before deciding to drop the rpy2/R backend in favor of pylimma.
Runs both backends on the same synthetic data (with a known injected
effect) and reports how closely they agree.

Requires BOTH backends installed: R + limma (for rpy2) AND pylimma.
See README.md for R/rpy2 setup.

Usage:
    python scripts/compare_limma_backends.py
"""
import numpy as np
import pandas as pd

from omics_app.stats.limma_bridge import run_limma_two_group
from omics_app.stats.pylimma_bridge import run_pylimma_two_group

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n_proteins = 200
    test_cols = ["test_1", "test_2", "test_3"]
    control_cols = ["control_1", "control_2", "control_3"]

    data = pd.DataFrame(
        rng.normal(10, 1, size=(n_proteins, 6)),
        columns=test_cols + control_cols,
        index=[f"PROT_{i}" for i in range(n_proteins)],
    )
    data.loc[data.index[:20], test_cols] += 3  # inject real fold change

    group_labels = ["test"] * 3 + ["control"] * 3

    r_result = run_limma_two_group(data, "test", "control", group_labels)
    py_result = run_pylimma_two_group(data, "test", "control", group_labels)

    r_result = r_result.add_suffix("_rpy2")
    py_result = py_result.add_suffix("_pylimma")
    combined = r_result.join(py_result)

    for col in ["logFC", "P.Value", "adj.P.Val"]:
        diff = (combined[f"{col}_rpy2"] - combined[f"{col}_pylimma"]).abs()
        print(f"{col:12s}  max abs diff: {diff.max():.2e}   mean abs diff: {diff.mean():.2e}")

    r_sig = set(combined.index[combined["P.Value_rpy2"] < 0.05])
    py_sig = set(combined.index[combined["P.Value_pylimma"] < 0.05])
    print(f"\nSignificant (p<0.05) -- rpy2: {len(r_sig)}, pylimma: {len(py_sig)}, "
          f"agree on: {len(r_sig & py_sig)}, disagree on: {len(r_sig ^ py_sig)}")

    if r_sig ^ py_sig:
        print("Proteins where the two backends disagree on significance:",
              sorted(r_sig ^ py_sig))
