# Multi-Omics Analysis Suite — Python Port

Python/Dash port of `app_12-02.R`, a Shiny app for proteomics differential expression +
pathway enrichment. See `data/sample/synthetic_proteomics.xlsx` (generate it with
`scripts/generate_sample_dataset.py`) for a synthetic dataset with known ground truth, used
throughout development to validate the pipeline against a known-correct answer.

## Project status

| Tab | Status |
|---|---|
| 🏠 Home | ✅ Built |
| 📁 Data Upload | ✅ Built |
| ⚖️ Comparisons | ✅ Built (2-group path only — see Known Gaps) |
| 🎨 Color Mapping | ✅ Built |
| 🔬 Analysis | ✅ Built (2-group path only — see Known Gaps) |
| 📊 Visualization | ❌ Not started |
| 🧬 Enrichment | ❌ Not started |
| 💾 Download | ❌ Not started |
| ℹ️ Session Info | ❌ Not started |

## Known gaps — read before trusting output on real data

1. **MBQN normalization needs more validation.** `stats/mbqn.py` implements the actual
   algorithm (row-median-center → quantile-normalize via `pylimma.normalize.normalize_quantiles`
   → add median back), replacing the earlier simple-median-centering placeholder. On the
   synthetic benchmark, swapping in real MBQN caused ground-truth recovery to collapse (35/40 →
   4/40 true positives, 0 → 18 false positives) compared to the placeholder. The cause isn't
   confirmed: `normalize_quantiles` itself checks out correctly in isolation (verified with a
   controlled test), and the algorithm structure matches MBQN's documented approach, but real
   R's `MBQN` package isn't installable in every environment for a direct side-by-side diff
   (unlike `limma`, it's not packaged for Ubuntu/apt, and CRAN/Bioconductor may not be reachable
   depending on your network setup). **Test this against your own real data and compare to the
   R app's output before trusting Analysis tab numbers.** If it doesn't hold up, reverting to
   simple median-centering (see git history for `_median_center_normalize` in `pipeline.py`) is
   a safer interim default.
2. **ANOVA (multi-group) comparisons aren't wired to a working pipeline.** Comparisons tab
   collects ANOVA group selections, but Analysis tab shows "not yet supported" for them — the
   ANOVA path doesn't have per-group abundance columns auto-derived from group flags the way
   the 2-group path does. `stats/anova_path.py` (pure scipy/statsmodels, no R needed) has the
   underlying stats logic; it just isn't connected to Comparisons/Analysis yet.
3. **No manual group-entry fallback.** Comparisons and Color Mapping both assume the uploaded
   file has "Found in Sample Group" columns. Data without those columns (R's manual-entry path)
   isn't supported yet.
4. **"Run All Comparisons" button has no callback.** Only single-comparison "Run Analysis" works.

## Validated and safe to trust

- **The core DEA statistics engine (`stats/pylimma_bridge.py`) has been directly validated
  against real R limma.** Installed R + the actual `limma` package, ran both on identical data,
  and diffed every value: **zero difference** in logFC, P.Value, and adj.P.Val (to 8 decimal
  places), and identical significance calls on all 30 test cases. See "R validation" below for
  how to reproduce this.
- Filtering and imputation logic (`data/filtering.py`) is unit-tested and matches R's formulas
  exactly (including the `max(1, floor(n * ratio))` filter threshold and `set.seed(1)`-equivalent
  fixed random seed for reproducible imputation).
- Column auto-detection (`data/columns.py`) is unit-tested against R's `build_main_data_index()`
  logic.

## 1. Setup

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .                   # makes `omics_app` importable from anywhere
```

The `pip install -e .` step matters — without it, `python src/omics_app/app.py` fails with
`ModuleNotFoundError: No module named 'omics_app'`, since Python only adds the script's own
folder to the import path, not `src/`.

In VS Code: **Cmd/Ctrl+Shift+P → "Python: Select Interpreter" → `.venv`** (should happen
automatically via `.vscode/settings.json`, but confirm in the status bar).

## 2. Run the app

```bash
python src/omics_app/app.py
```
Opens on `http://127.0.0.1:8050`. Generate the synthetic test dataset first if you don't have
real data yet:
```bash
python scripts/generate_sample_dataset.py
```

## 3. Run the tests

```bash
pytest
```

## 4. The stats backend: pylimma (primary, no R needed)

`stats/pylimma_bridge.py` is the active DEA engine — a pure-Python limma port, validated
against real R limma (see above). No R installation is required for normal development or
running the app.

`stats/limma_bridge.py` (an `rpy2` → real R `limma` bridge) still exists as optional reference
tooling, **not used by the active pipeline** (`pipeline.py` only imports `pylimma_bridge`). Keep
it around if you want to re-run a direct comparison against real R limma on your own data:

```bash
python scripts/compare_limma_backends.py
```

This requires R + `limma` + `rpy2` installed (see below) — only needed if you want to
re-validate, not for normal use.

### Installing R + limma + rpy2 (optional, only for re-validation)

1. Install R: macOS `brew install r`, Ubuntu/Debian `sudo apt install r-base-core`, Windows from
   CRAN.
2. Install `limma`. On Debian/Ubuntu, a prebuilt package is available and avoids needing CRAN/
   Bioconductor network access:
   ```bash
   sudo apt install r-bioc-limma
   ```
   Otherwise, from an R console: `BiocManager::install("limma")`.
3. Install `rpy2`: `pip install rpy2`. Note: in some environments `rpy2`'s compiled bindings can
   hit a native ABI mismatch against the installed R build (`undefined symbol: R_getVar` or
   similar). If that happens, the most reliable workaround is bypassing `rpy2` entirely — call
   `Rscript` as a subprocess and exchange data via CSV, which is what `validate_limma_bridge.py`
   falls back to conceptually if you adapt it that way.
4. Run `python scripts/validate_limma_bridge.py` to confirm the bridge works before relying on
   it.

## Project layout

```
omics_port/
├── .vscode/
├── src/omics_app/
│   ├── app.py              # Dash entry point, sidebar-nav shell, all tabs mounted permanently
│   ├── server.py           # Flask-Caching instance (large file uploads, not client-side stores)
│   ├── data/
│   │   ├── filtering.py    # filter_valids, impute_downshift -- unit tested
│   │   └── columns.py      # build_main_data_index, extract_group_names_from_columns
│   ├── stats/
│   │   ├── pipeline.py         # orchestrates the full 2-group DEA pipeline (active path)
│   │   ├── pylimma_bridge.py   # pure-Python limma -- validated against real R limma
│   │   ├── mbqn.py             # real MBQN normalization -- needs more validation, see Known Gaps
│   │   ├── anova_path.py       # pure scipy/statsmodels ANOVA -- not yet wired to UI
│   │   └── limma_bridge.py     # rpy2 -> real R limma -- optional reference/re-validation only
│   ├── ui/
│   │   ├── home.py, upload.py, comparisons.py, colors.py, analysis.py
│   ├── enrichment/, plotting/, export/    # empty stubs, not started
├── scripts/
│   ├── generate_sample_dataset.py     # synthetic data with known ground truth
│   ├── validate_limma_bridge.py       # sanity-checks the rpy2 bridge (optional)
│   └── compare_limma_backends.py      # diffs rpy2 vs pylimma on synthetic data (optional)
├── tests/
└── requirements.txt
```

## State management note

The R app's `rv` (`reactiveValues`) maps to Dash `dcc.Store` components in `app.py`. Uploaded
file bytes go through a server-side cache (`server.py`, `flask-caching`) rather than a
client-side store, since the full base64 file content round-tripping through the browser on
every callback was the actual cause of slow uploads early on — only a small cache token lives
in the client-side store now.

All tab panels are built once at startup and stay permanently in the DOM; switching tabs only
toggles `display: block`/`none` (see `set_active_panel` in `app.py`). This is deliberate —
destroying and rebuilding tab content on every switch was causing dropdown selections and
loaded data to reset when navigating away and back.