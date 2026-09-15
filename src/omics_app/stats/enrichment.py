"""
Port of run_enrichment_analysis (app_12-02.R, ~lines 9700-9950) --
Bioconductor engine only (clusterProfiler::enrichGO +
ReactomePA::enrichPathway, run separately on up- and down-regulated
significant proteins). The DAVID Webservice engine (R's other
enrichment_engine option) isn't ported; ui/enrichment.py's DAVID
settings panel stays disabled and the run callback shows a message
instead of attempting a DAVID call.

R's engine needs Bioconductor organism-annotation packages
(org.Hs.eg.db, org.Mm.eg.db, ...) with no direct Python equivalent, so
this uses Enrichr (via gseapy) instead: a web-API-based
over-representation test against curated GO/Reactome gene sets.
Practically this changes two things from R's approach:
  - No Entrez ID mapping step (R's bitr() call, lines 9713/9749) --
    Enrichr works directly off gene SYMBOLS, so R's id-mapping logic
    just doesn't apply here.
  - Requires network access at analysis time (calls
    maayanlab.cloud's Enrichr API through gseapy) -- R's local
    Bioconductor packages don't need network once installed. If the
    app is deployed somewhere without outbound internet, this will
    fail with no offline fallback.

Organism coverage: Enrichr only has dedicated instances for Human,
Mouse, Fly, Yeast, Worm, and Fish (the "modEnrichr" federation). Rat
and Arabidopsis (2 of ui/enrichment.py's 8 ORGANISM_OPTIONS) have no
Enrichr instance, so they fall back to the Human gene-set library --
results for those two should be treated as approximate, and the
returned "warning" field says so.

Dependency: `gseapy` (add to requirements.txt).
"""

import pandas as pd

try:
    import gseapy as gp
except ImportError:  # pragma: no cover - exercised only when gseapy isn't installed
    gp = None

GENE_SETS = {
    "GO_Biological_Process": "GO_Biological_Process_2023",
    "GO_Cellular_Component": "GO_Cellular_Component_2023",
    "GO_Molecular_Function": "GO_Molecular_Function_2023",
    "Reactome_Pathways": "Reactome_2022",
}

# Enrichr/modEnrichr only has dedicated instances for these organisms.
_ORGANISM_TO_ENRICHR = {
    "human": "human",
    "mouse": "mouse",
    "fly": "fly",
    "worm": "worm",
    "yeast": "yeast",
    "zebrafish": "fish",
    "rat": "human",
    "arabidopsis": "human",  # no native instance -- approximated
}
_NO_NATIVE_SUPPORT = {"rat", "arabidopsis"}


class EnrichmentError(Exception):
    """User-facing validation/failure, caught by the UI layer -- mirrors
    AnalysisError in stats/pipeline.py."""


def _gene_symbols_for_rows(df: pd.DataFrame) -> pd.Series:
    """Minimal, enrichment-specific gene-symbol resolution: prefers an
    explicit Gene Symbol column, else parses GN=... out of Description."""
    if "Gene Symbol" in df.columns:
        symbol = df["Gene Symbol"].astype(str)
    elif "Gene.Symbol" in df.columns:
        symbol = df["Gene.Symbol"].astype(str)
    elif "Description" in df.columns:
        symbol = df["Description"].astype(str).str.extract(r" GN=(\S+)")[0]
    else:
        symbol = pd.Series(pd.NA, index=df.index)
    return symbol.astype(str).str.strip()


def _gene_list(df: pd.DataFrame, sig_mask: pd.Series, direction: str) -> list[str]:
    """direction: 'up' or 'down', split by logFC sign among significant rows."""
    sig = df[sig_mask.reindex(df.index, fill_value=False)]
    if sig.empty or "logFC" not in sig.columns:
        return []
    subset = sig[sig["logFC"] > 0] if direction == "up" else sig[sig["logFC"] < 0]
    if subset.empty:
        return []
    symbols = _gene_symbols_for_rows(subset)
    symbols = symbols[symbols.notna() & (symbols != "") & (symbols.str.lower() != "nan")]
    return sorted(set(symbols))


def run_enrichment_analysis(
    df: pd.DataFrame,
    sig_mask: pd.Series,
    organism: str,
    gene_set_keys: list[str],
    pvalue_cutoff: float = 0.05,
    qvalue_cutoff: float = 1.0,
) -> dict:
    """
    df, sig_mask: same shape as stats/comparison_data.py's
                  load_comparison() output (df = comparison results,
                  sig_mask = boolean Series over df.index for the
                  chosen significance column).
    organism: one of ui/enrichment.py's ORGANISM_OPTIONS values.
    gene_set_keys: subset of GENE_SETS keys to run (ui/enrichment.py's
                  "GO" checklist value expands to all 3 GO aspects,
                  "Reactome" to Reactome_Pathways).

    Returns {"results_df", "up_genes", "down_genes", "warning"}.
    """
    if gp is None:
        raise EnrichmentError(
            "The 'gseapy' package isn't installed -- run `pip install gseapy` first."
        )

    up_genes = _gene_list(df, sig_mask, "up")
    down_genes = _gene_list(df, sig_mask, "down")
    if not up_genes and not down_genes:
        raise EnrichmentError(
            "No significant proteins with resolvable gene symbols for this comparison."
        )

    enrichr_organism = _ORGANISM_TO_ENRICHR.get(organism, "human")
    warning = None
    if organism in _NO_NATIVE_SUPPORT:
        warning = (
            f"Enrichr has no dedicated instance for this organism; results were computed "
            f"against the {enrichr_organism} gene sets and should be treated as approximate."
        )

    gene_sets = [GENE_SETS[k] for k in gene_set_keys if k in GENE_SETS] or list(GENE_SETS.values())

    frames = []
    for direction, genes in [("Up", up_genes), ("Down", down_genes)]:
        if not genes:
            continue
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=gene_sets,
                organism=enrichr_organism,
                outdir=None,
                no_plot=True,
            )
        except Exception as exc:
            raise EnrichmentError(
                f"Enrichr request failed for {direction}-regulated genes: {exc}"
            ) from exc

        result = enr.results.copy()
        result["Direction"] = direction
        frames.append(result)

    if not frames:
        raise EnrichmentError("No enrichment results returned.")

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.rename(
        columns={
            "Gene_set": "Source",
            "P-value": "PValue",
            "Adjusted P-value": "AdjPValue",
        }
    )
    keep_cols = [
        c
        for c in [
            "Direction",
            "Source",
            "Term",
            "Overlap",
            "PValue",
            "AdjPValue",
            "Odds Ratio",
            "Genes",
        ]
        if c in combined.columns
    ]
    combined = combined[keep_cols]
    combined = (
        combined[(combined["PValue"] <= pvalue_cutoff) & (combined["AdjPValue"] <= qvalue_cutoff)]
        .sort_values("PValue")
        .reset_index(drop=True)
    )

    return {
        "results_df": combined,
        "up_genes": up_genes,
        "down_genes": down_genes,
        "warning": warning,
    }
