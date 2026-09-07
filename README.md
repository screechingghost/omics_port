# Setup

## 1. Open the project correctly
Open the **`omics_port/` folder itself** in VS Code (File → Open Folder), not a parent folder —
the `.vscode/` settings, and the `src` layout (`python.analysis.extraPaths`) assume the folder
root is the workspace root.

## 2. Install the recommended extensions
VS Code will prompt you automatically ("This workspace has extension recommendations") because
of `.vscode/extensions.json`. Accept it, or open the Extensions panel and search
`@recommended`. Key ones:
- **Python** + **Pylance** + **debugpy** — core Python support
- **Ruff** — linting/formatting (configured to run on save)
- **R** (REditorSupport.r) — only needed if you want to inspect/run R snippets directly in VS
  Code while porting; not required to run the app itself.
- **Jupyter** — handy for interactively comparing R vs. Python output row-by-row during
  validation.

## 3. Create the virtual environment
```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .                   # makes `omics_app` importable from anywhere
```
The `pip install -e .` step is important — without it, running `python src/omics_app/app.py`
directly fails with `ModuleNotFoundError: No module named 'omics_app'`, because Python only adds
the script's own folder to the import path, not `src/`. The editable install registers the
package properly so it works no matter how or where you run it.

Then in VS Code: **Cmd/Ctrl+Shift+P → "Python: Select Interpreter" → `.venv`**. This should
happen automatically given `python.defaultInterpreterPath` in settings, but confirm it in the
bottom-right status bar — a wrong interpreter is the single most common source of "why can't it
find my package" confusion.

## 4. R + rpy2 setup (only needed for the limma / 2-group path)
This is the one non-obvious part. `rpy2` needs an actual R installation on your machine —
installing the `rpy2` pip package alone is not enough.

1. Install R (if not already installed):
   - macOS: `brew install r`
   - Ubuntu/Debian: `sudo apt install r-base`
   - Windows: download from https://cran.r-project.org/bin/windows/base/
2. Install limma inside R itself (one-time, from an R console — not pip):
   ```r
   install.packages("BiocManager")
   BiocManager::install("limma")
   ```
3. Confirm `rpy2` can find your R installation. `.vscode/settings.json` sets `R_HOME` for
   Linux/macOS — **edit the path if yours differs** (check with `R RHOME` in a terminal).
   On Windows, add `R_HOME` and put `R.dll`'s folder on `PATH` instead; rpy2's own install docs
   cover the Windows-specific quirks in more detail than is worth duplicating here.
4. Run the validation script before writing anything else that depends on this bridge:
   ```bash
   python scripts/validate_limma_bridge.py
   ```
   or use the **"Validate: limma bridge (rpy2)"** debug config (F5, pick it from the dropdown).
   It runs limma on synthetic data with a known injected effect and reports how many of the
   20 "real" signals came back significant — if this doesn't run cleanly, fix it now, because
   every later step in the 2-group path depends on it.

**Note:** the ANOVA (multi-group) path does *not* need any of this — it's pure
`scipy`/`statsmodels`. Only skip this section entirely if you're certain you won't touch the
2-group comparison path yet.

## 4b. Optional: pylimma instead of rpy2 (validate first)
`pylimma` (`stats/pylimma_bridge.py`) is a pure-Python limma port -- no R install needed at all,
which would remove the entire rpy2/R setup above. It's a young package (first release May 2026),
so validate it against the R backend on your real data before relying on it alone:
```bash
python scripts/compare_limma_backends.py
```
This requires *both* backends installed (R+limma+rpy2, and `pip install pylimma`) since it runs
the same synthetic data through each and reports the max/mean difference in logFC, P.Value, and
adj.P.Val, plus which proteins (if any) the two backends disagree on for significance. If the
differences are within a tolerance you're comfortable with on your actual dataset, switch to
`pylimma_bridge.py` and drop the R/rpy2 setup entirely; if not, keep the rpy2 path as your
source of truth.

Note pylimma's `top_table()` returns snake_case columns (`log_fc`, `p_value`, `adj_p_value`),
already renamed to R's convention (`logFC`, `P.Value`, `adj.P.Val`) inside `pylimma_bridge.py`
so the rest of the codebase doesn't need to know which backend is in use.

## 5. Run the tests
```bash
pytest
```
or use the **"Pytest: Current file"** debug config with a test file open, or the Testing panel
(flask icon) in the sidebar, which Pylance/pytest integration populates automatically once the
interpreter is set correctly.

## 6. Run the app shell
Press **F5** with **"Dash: Run app"** selected, or:
```bash
python src/omics_app/app.py
```
Opens on `http://127.0.0.1:8050` with the 9 tabs stubbed out. Fill in `render_tab()` in
`src/omics_app/app.py` as you port each tab.

## Project layout
```
omics_port/
├── .vscode/              # settings, debug configs, extension recommendations
├── src/omics_app/
│   ├── app.py            # Dash entry point, 9-tab shell
│   ├── data/             # filtering, imputation, column parsing (port first)
│   ├── stats/
│   │   ├── anova_path.py     # pure scipy/statsmodels, no R needed
│   │   └── limma_bridge.py   # rpy2 -> R limma, 2-group path only
│   ├── enrichment/        # gseapy + DAVID client go here
│   ├── plotting/          # plotly/seaborn plot generators
│   └── export/            # Excel/PDF/zip export
├── scripts/
│   └── validate_limma_bridge.py
├── tests/
└── requirements.txt
```

## State management note
The R app's `rv` (`reactiveValues`) becomes Dash `dcc.Store` components in `app.py`. Dash's
default store serializes to JSON on the client side — fine for small config values (comparison
definitions, color mapping) but likely too slow/large for the actual data matrix and DEA results
if your uploads approach the 200MB limit the R app allows. If you hit that wall, switch those
specific stores to `storage_type='memory'` backed by `flask-caching` on the server side rather
than trying to push large DataFrames through the browser.
