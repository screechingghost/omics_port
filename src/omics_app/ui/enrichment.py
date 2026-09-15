"""
Port of R tabPanel("🧬 Enrichment", ...) (app_12-02.R, lines 1385-1502)
and observeEvent(input$run_enrichment_btn, ...) (~line 9700) --
Bioconductor engine only. See stats/enrichment.py for the actual
Enrichr-based over-representation logic (UI-free, testable) and
plotting/enrichment.py for the 2 wired plot types; this module is the
thin Dash wrapper around both.

What's wired:
  - Comparison selector, driven by store-dea-results (only comparisons
    with a completed 2-group Analysis run show up).
  - Run Enrichment Analysis -> stats/enrichment.run_enrichment_analysis
    -> store-enrichment-results, keyed by comparison name (same
    accumulate-don't-overwrite pattern ui/analysis.py uses for
    store-dea-results).
  - Results summary, results table (dash_table.DataTable), and 2 of
    the 9 plot types (Grouped Barplot, Dotplot).
  - Analysis Log entry on every successful run (store-analysis-log).

What's still NOT wired:
  - DAVID Webservice engine -- selecting it and clicking Run just
    shows a message; none of the david-* inputs are read anywhere.
  - "Run All Enrichment Analyses" button -- stays disabled.
  - 7 of 9 plot types -- selecting one shows a placeholder message.
  - "Download TIFF" for the enrichment plot -- stays disabled.
  - "All Plot Preview" section at the bottom -- stays a static
    placeholder.
"""

# pyright: reportCallIssue=false, reportInvalidTypeForm=false
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, callback, dash_table, dcc, html, no_update
from dash.exceptions import PreventUpdate

from omics_app.plotting.enrichment import build_enrichment_barplot, build_enrichment_dotplot
from omics_app.stats.comparison_data import load_comparison
from omics_app.stats.enrichment import EnrichmentError, run_enrichment_analysis
from omics_app.stats.session_log import append_log_entry

ORGANISM_OPTIONS = [
    {"label": "Human (GO + Reactome)", "value": "human"},
    {"label": "Mouse (GO + Reactome)", "value": "mouse"},
    {"label": "Rat (approximated via Human gene sets)", "value": "rat"},
    {"label": "Zebrafish (GO + Reactome)", "value": "zebrafish"},
    {"label": "Fly (GO only)", "value": "fly"},
    {"label": "Worm (C. elegans) (GO only)", "value": "worm"},
    {"label": "Yeast (GO + Reactome)", "value": "yeast"},
    {"label": "Arabidopsis (approximated via Human gene sets)", "value": "arabidopsis"},
]

PLOT_TYPE_OPTIONS = [
    {"label": "Grouped Barplot", "value": "barplot"},
    {"label": "Dotplot", "value": "dotplot"},
    {"label": "GO: By Ontology (BP/MF/CC)", "value": "go_ontology"},
    {"label": "Bubble Plot", "value": "bubble"},
    {"label": "Up vs Down Comparison", "value": "up_down_compare"},
    {"label": "Gene-Pathway Network (cnetplot)", "value": "cnetplot"},
    {"label": "Enrichment Map (term similarity)", "value": "emapplot"},
    {"label": "Chord: Pathway Up", "value": "chord_up"},
    {"label": "Chord: Pathway Down", "value": "chord_down"},
]

_WIRED_PLOT_TYPES = {"barplot", "dotplot"}


def _empty_figure(message: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font=dict(size=14))
    fig.update_layout(xaxis_visible=False, yaxis_visible=False, template="plotly_white")
    return fig


def layout() -> html.Div:
    settings = dbc.Card(
        dbc.CardBody(
            [
                html.H4("Enrichment Analysis Setup"),
                dbc.Label("Comparison:", className="mt-2"),
                dcc.Dropdown(
                    id="enrichment-comparison-select", options=[], placeholder="Run Analysis first"
                ),
                dbc.Label("Select Organism:", className="mt-3"),
                dcc.Dropdown(
                    id="enrichment-organism",
                    options=ORGANISM_OPTIONS,
                    value="human",
                    clearable=False,
                ),
                dbc.Label("Enrichment Engine:", className="mt-3"),
                dbc.RadioItems(
                    id="enrichment-engine",
                    options=[
                        {"label": "Bioconductor GO/Reactome (via Enrichr)", "value": "bioc"},
                        {"label": "DAVID Webservice (not wired up)", "value": "david"},
                    ],
                    value="bioc",
                ),
                html.Div(
                    [
                        html.Hr(),
                        html.H6("Bioconductor settings", style={"color": "var(--muted)"}),
                        dbc.Label("Analysis Types:"),
                        dbc.Checklist(
                            id="enrichment-types",
                            options=[
                                {"label": "Gene Ontology (GO)", "value": "GO"},
                                {"label": "Reactome Pathways", "value": "Reactome"},
                            ],
                            value=["GO", "Reactome"],
                        ),
                    ]
                ),
                html.Div(
                    [
                        html.Hr(),
                        html.H6("DAVID settings (unused)", style={"color": "var(--muted)"}),
                        dbc.Label("DAVID Email:"),
                        dbc.Input(id="david-email", type="email", placeholder="you@lab.edu"),
                        dbc.Label("DAVID Species:", className="mt-2"),
                        dbc.Input(id="david-species", type="text", value="Homo sapiens"),
                        dbc.Label("DAVID ID Type:", className="mt-2"),
                        dcc.Dropdown(
                            id="david-id-type",
                            options=[
                                {"label": "Gene Symbol", "value": "OFFICIAL_GENE_SYMBOL"},
                                {"label": "Entrez Gene ID", "value": "ENTREZ_GENE_ID"},
                                {"label": "UniProt Accession", "value": "UNIPROT_ACCESSION"},
                                {"label": "RefSeq Protein", "value": "REFSEQ_PROTEIN_ACCESSION"},
                            ],
                            value="OFFICIAL_GENE_SYMBOL",
                            clearable=False,
                        ),
                        dbc.Label("DAVID Categories:", className="mt-2"),
                        dbc.Input(
                            id="david-categories",
                            type="text",
                            value="GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,KEGG_PATHWAY",
                        ),
                        dbc.Label("DAVID EASE threshold:", className="mt-2"),
                        dbc.Input(
                            id="david-ease-threshold",
                            type="number",
                            value=1,
                            min=0.001,
                            max=1,
                            step=0.005,
                        ),
                        dbc.Label("DAVID minimum count:", className="mt-2"),
                        dbc.Input(
                            id="david-count-threshold",
                            type="number",
                            value=1,
                            min=1,
                            max=100,
                            step=1,
                        ),
                        dbc.Checklist(
                            id="david-mock",
                            options=[
                                {"label": "Use mock DAVID for local testing only", "value": "mock"}
                            ],
                            value=[],
                            switch=True,
                            className="mt-2",
                        ),
                    ]
                ),
                html.Hr(),
                dbc.Label("Enrichment p-value cutoff:"),
                dbc.Input(
                    id="enrichment-pvalue-cutoff",
                    type="number",
                    value=0.05,
                    min=0.001,
                    max=1,
                    step=0.005,
                ),
                dbc.Label("Enrichment q-value cutoff:", className="mt-2"),
                dbc.Input(
                    id="enrichment-qvalue-cutoff",
                    type="number",
                    value=1,
                    min=0.001,
                    max=1,
                    step=0.005,
                ),
                dbc.Button(
                    "🧬 Run Enrichment Analysis",
                    id="run-enrichment-btn",
                    color="primary",
                    size="lg",
                    className="w-100 mt-3",
                    disabled=True,
                ),
                dbc.Button(
                    "▶ Run All Enrichment Analyses",
                    id="run-all-enrichment-btn",
                    color="info",
                    className="w-100 mt-2",
                    disabled=True,
                ),
                html.Hr(),
                html.H5("Progress:"),
                html.Pre(
                    "Not run yet.",
                    id="enrichment-progress",
                    style={
                        "fontSize": "12px",
                        "color": "var(--muted)",
                        "backgroundColor": "var(--bg)",
                        "padding": "8px",
                        "borderRadius": "6px",
                        "whiteSpace": "pre-wrap",
                    },
                ),
            ]
        ),
        className="card-accent-teal",
    )

    results = dbc.Card(
        dbc.CardBody(
            [
                html.H4("Enrichment Results"),
                html.Pre(
                    "Run an enrichment analysis to see a summary here.",
                    id="enrichment-summary",
                    style={"fontSize": "13px", "color": "var(--muted)", "whiteSpace": "pre-wrap"},
                ),
                html.Hr(),
                html.H5("Enrichment Table:"),
                html.Div(
                    "No results yet.",
                    id="enrichment-results-table",
                    style={"color": "var(--muted)", "fontSize": "13px"},
                ),
                html.Hr(),
                dbc.Row(
                    [
                        dbc.Col(html.H5("Enrichment Visualization:"), width=6),
                        dbc.Col(
                            dcc.Dropdown(
                                id="enrichment-plot-type",
                                options=PLOT_TYPE_OPTIONS,
                                value="barplot",
                                clearable=False,
                            ),
                            width=6,
                        ),
                    ]
                ),
                dcc.Loading(
                    dcc.Graph(
                        id="plot-enrichment",
                        figure=_empty_figure("Run an enrichment analysis first."),
                        style={"height": "440px", "marginTop": "12px"},
                    )
                ),
                dbc.Button(
                    "⬇ Download TIFF",
                    id="download-enrichment-plot-btn",
                    color="success",
                    size="sm",
                    outline=True,
                    className="mt-2",
                    disabled=True,
                ),
                html.Hr(),
                html.H5("All Plot Preview"),
                html.Div(
                    "Not wired up yet.",
                    id="enrichment-all-plots-preview",
                    style={"color": "var(--muted)", "fontSize": "13px"},
                ),
            ]
        ),
        className="card-accent-indigo",
    )

    return dbc.Row([dbc.Col(settings, width=4), dbc.Col(results, width=8)])


@callback(
    Output("enrichment-comparison-select", "options"),
    Output("enrichment-comparison-select", "value"),
    Input("store-dea-results", "data"),
)
def populate_enrichment_comparison_selector(dea_results):
    """Mirrors ui/visualization.py's comparison-selector callback, but
    also filters out comparisons load_comparison() can't reconstruct
    (ANOVA / pre-upgrade entries)."""
    if not dea_results:
        return [], None
    valid_names = [name for name in dea_results if load_comparison(dea_results, name) is not None]
    if not valid_names:
        return [], None
    return [{"label": name, "value": name} for name in valid_names], valid_names[0]


@callback(
    Output("run-enrichment-btn", "disabled"),
    Input("enrichment-comparison-select", "value"),
)
def toggle_run_enrichment_button(comp_name):
    return not comp_name


@callback(
    Output("store-enrichment-results", "data"),
    Output("enrichment-progress", "children"),
    Output("store-analysis-log", "data", allow_duplicate=True),
    Input("run-enrichment-btn", "n_clicks"),
    State("enrichment-comparison-select", "value"),
    State("enrichment-engine", "value"),
    State("enrichment-organism", "value"),
    State("enrichment-types", "value"),
    State("enrichment-pvalue-cutoff", "value"),
    State("enrichment-qvalue-cutoff", "value"),
    State("store-dea-results", "data"),
    State("store-enrichment-results", "data"),
    State("store-analysis-log", "data"),
    prevent_initial_call=True,
)
def run_enrichment(
    n_clicks,
    comp_name,
    engine,
    organism,
    types,
    pval_cutoff,
    qval_cutoff,
    dea_results,
    existing_enrichment_results,
    existing_log,
):
    """Port of observeEvent(input$run_enrichment_btn, ...) (R ~line 9700),
    Bioconductor engine only."""
    if engine != "bioc":
        return (
            no_update,
            "DAVID Webservice engine isn't wired up yet -- switch to Bioconductor GO/Reactome to run enrichment.",
            no_update,
        )

    loaded = load_comparison(dea_results, comp_name)
    if loaded is None:
        return (
            no_update,
            "Select a completed 2-group comparison first (run Analysis, then come back here).",
            no_update,
        )

    types = types or []
    gene_set_keys = []
    if "GO" in types:
        gene_set_keys += ["GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function"]
    if "Reactome" in types:
        gene_set_keys.append("Reactome_Pathways")
    if not gene_set_keys:
        return (
            no_update,
            "Select at least one of GO or Reactome under Bioconductor settings.",
            no_update,
        )

    try:
        result = run_enrichment_analysis(
            loaded["df"],
            loaded["sig_mask"],
            organism or "human",
            gene_set_keys,
            pvalue_cutoff=pval_cutoff or 0.05,
            qvalue_cutoff=qval_cutoff or 1.0,
        )
    except EnrichmentError as exc:
        return no_update, f"Enrichment failed: {exc}", no_update

    lines = [
        f"Comparison: {comp_name}",
        f"Up-regulated genes tested: {len(result['up_genes'])}",
        f"Down-regulated genes tested: {len(result['down_genes'])}",
        f"Enriched terms passing cutoffs: {len(result['results_df'])}",
    ]
    if result["warning"]:
        lines.append(f"Note: {result['warning']}")

    payload = dict(existing_enrichment_results or {})
    payload[comp_name] = {
        "results": result["results_df"].to_dict("records"),
        "up_genes": result["up_genes"],
        "down_genes": result["down_genes"],
        "warning": result["warning"],
    }
    log = append_log_entry(
        existing_log,
        f"Enrichment run: '{comp_name}' -- {len(result['results_df'])} term(s) passing cutoffs",
    )
    return payload, "\n".join(lines), log


@callback(
    Output("enrichment-summary", "children"),
    Input("store-enrichment-results", "data"),
    Input("enrichment-comparison-select", "value"),
)
def render_enrichment_summary(enrichment_results, comp_name):
    if not enrichment_results or not comp_name or comp_name not in enrichment_results:
        raise PreventUpdate
    entry = enrichment_results[comp_name]
    lines = [
        f"Comparison: {comp_name}",
        f"Up-regulated genes tested: {len(entry.get('up_genes', []))}",
        f"Down-regulated genes tested: {len(entry.get('down_genes', []))}",
        f"Enriched terms: {len(entry.get('results', []))}",
    ]
    if entry.get("warning"):
        lines.append(f"Note: {entry['warning']}")
    return "\n".join(lines)


@callback(
    Output("enrichment-results-table", "children"),
    Input("store-enrichment-results", "data"),
    Input("enrichment-comparison-select", "value"),
)
def render_enrichment_table(enrichment_results, comp_name):
    if not enrichment_results or not comp_name or comp_name not in enrichment_results:
        raise PreventUpdate
    records = enrichment_results[comp_name].get("results", [])
    if not records:
        return "No enriched terms passed the cutoffs."

    df = pd.DataFrame(records)
    display_cols = [
        c
        for c in ["Direction", "Source", "Term", "Overlap", "PValue", "AdjPValue", "Genes"]
        if c in df.columns
    ]
    df = df[display_cols].head(200).copy()
    for col in ("PValue", "AdjPValue"):
        if col in df.columns:
            df[col] = df[col].map(lambda v: f"{v:.2e}")
    if "Genes" in df.columns:
        df["Genes"] = df["Genes"].astype(str).str.slice(0, 80)

    return dash_table.DataTable(
        columns=[{"name": c, "id": c} for c in display_cols],
        data=df.to_dict("records"),
        page_size=15,
        sort_action="native",
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_cell={
            "fontSize": "12px",
            "fontFamily": "Inter, sans-serif",
            "textAlign": "left",
            "padding": "6px",
        },
        style_header={"fontWeight": "600", "backgroundColor": "var(--bg)"},
    )


@callback(
    Output("plot-enrichment", "figure"),
    Input("store-enrichment-results", "data"),
    Input("enrichment-comparison-select", "value"),
    Input("enrichment-plot-type", "value"),
)
def render_enrichment_plot(enrichment_results, comp_name, plot_type):
    if not enrichment_results or not comp_name or comp_name not in enrichment_results:
        raise PreventUpdate
    records = enrichment_results[comp_name].get("results", [])
    df = pd.DataFrame(records)
    if df.empty:
        return _empty_figure("No enriched terms passed the cutoffs.")

    if plot_type not in _WIRED_PLOT_TYPES:
        label = next((o["label"] for o in PLOT_TYPE_OPTIONS if o["value"] == plot_type), plot_type)
        return _empty_figure(f"'{label}' isn't wired up yet -- try Grouped Barplot or Dotplot.")

    if plot_type == "dotplot":
        return build_enrichment_dotplot(df)
    return build_enrichment_barplot(df)
