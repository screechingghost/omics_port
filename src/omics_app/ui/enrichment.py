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
  - Results summary, results table (dash_table.DataTable), and all 9
    plot types from plotting/enrichment.py (2 direct ports, 3 adapted to
    this port's flat Enrichr results shape, 4 built as approximations
    where R depends on packages/layouts with no wired Python equivalent
    -- see plotting/enrichment.py's module docstring for exactly which
    is which).
  - Download TIFF for the currently-shown enrichment plot.
  - "All Plot Preview" -- a 2-column grid rendering all 9 plot types at
    once for the selected comparison (a convenience addition; R has no
    literal equivalent of this specific grid).
  - Analysis Log entry on every successful run (store-analysis-log).

What's still NOT wired:
  - DAVID Webservice engine -- selecting it and clicking Run just
    shows a message; none of the david-* inputs are read anywhere.
  - "Run All Enrichment Analyses" button -- stays disabled.
"""

# pyright: reportCallIssue=false, reportInvalidTypeForm=false
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, callback, dash_table, dcc, html, no_update
from dash.exceptions import PreventUpdate

from omics_app.plotting.enrichment import (
    build_enrichment_barplot,
    build_enrichment_bubble_plot,
    build_enrichment_chord_plot,
    build_enrichment_cnetplot,
    build_enrichment_dotplot,
    build_enrichment_map_plot,
    build_go_ontology_plot,
    build_up_down_comparison_plot,
)
from omics_app.stats.comparison_data import load_comparison
from omics_app.stats.enrichment import EnrichmentError, run_enrichment_analysis
from omics_app.stats.session_log import append_log_entry
from omics_app.ui.visualization import _figure_to_tiff_bytes

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

_WIRED_PLOT_TYPES = {
    "barplot",
    "dotplot",
    "go_ontology",
    "bubble",
    "up_down_compare",
    "cnetplot",
    "emapplot",
    "chord_up",
    "chord_down",
}


def _empty_figure(message: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font={"size": 14})
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
                dcc.Download(id="download-enrichment-plot"),
                html.Hr(),
                html.H5("All Plot Preview"),
                dcc.Loading(
                    html.Div(
                        "Run an enrichment analysis to see all plot types here.",
                        id="enrichment-all-plots-preview",
                        style={"color": "var(--muted)", "fontSize": "13px"},
                    )
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

    return _build_plot_by_type(df, plot_type)


def _build_plot_by_type(df: pd.DataFrame, plot_type: str) -> go.Figure:
    """Dispatch table matching R's generate_enrichment_plot_by_type
    switch (R lines 10589-10598) -- all 9 plot_type options, see
    plotting/enrichment.py for what each one actually does (2 direct
    ports, 3 adapted to this port's flat results shape, 4 approximations
    where R relies on packages/layouts with no Python equivalent here)."""
    if plot_type == "dotplot":
        return build_enrichment_dotplot(df)
    if plot_type == "go_ontology":
        return build_go_ontology_plot(df)
    if plot_type == "bubble":
        return build_enrichment_bubble_plot(df)
    if plot_type == "up_down_compare":
        return build_up_down_comparison_plot(df)
    if plot_type == "cnetplot":
        return build_enrichment_cnetplot(df)
    if plot_type == "emapplot":
        return build_enrichment_map_plot(df)
    if plot_type == "chord_up":
        return build_enrichment_chord_plot(df, "up")
    if plot_type == "chord_down":
        return build_enrichment_chord_plot(df, "down")
    return build_enrichment_barplot(df)


@callback(
    Output("download-enrichment-plot-btn", "disabled"),
    Input("store-enrichment-results", "data"),
    Input("enrichment-comparison-select", "value"),
)
def toggle_download_enrichment_button(enrichment_results, comp_name):
    return not (enrichment_results and comp_name and comp_name in enrichment_results)


@callback(
    Output("download-enrichment-plot", "data"),
    Input("download-enrichment-plot-btn", "n_clicks"),
    State("plot-enrichment", "figure"),
    State("enrichment-comparison-select", "value"),
    State("enrichment-plot-type", "value"),
    prevent_initial_call=True,
)
def download_enrichment_tiff(n_clicks, figure, comp_name, plot_type):
    """Port of R's TIFF export for the enrichment plot, same pattern as
    ui/visualization.py's download_volcano_tiff etc."""
    if not figure:
        raise PreventUpdate
    tiff_bytes = _figure_to_tiff_bytes(figure)
    return dcc.send_bytes(
        lambda buf: buf.write(tiff_bytes),
        f"{comp_name or 'comparison'}_enrichment_{plot_type or 'plot'}.tiff",
    )


@callback(
    Output("enrichment-all-plots-preview", "children"),
    Input("store-enrichment-results", "data"),
    Input("enrichment-comparison-select", "value"),
)
def render_all_plots_preview(enrichment_results, comp_name):
    """Small-multiples grid of every wired plot type for the current
    comparison -- lets you compare all 9 views at once instead of
    switching the dropdown one at a time. R doesn't have a literal
    equivalent of this grid; it's a convenience addition on top of the
    ported single-plot dropdown above."""
    if not enrichment_results or not comp_name or comp_name not in enrichment_results:
        return "Run an enrichment analysis to see all plot types here."
    records = enrichment_results[comp_name].get("results", [])
    df = pd.DataFrame(records)
    if df.empty:
        return "No enriched terms passed the cutoffs."

    cards = []
    for opt in PLOT_TYPE_OPTIONS:
        plot_type = opt["value"]
        fig = _build_plot_by_type(df, plot_type)
        fig.update_layout(height=280, margin={"l": 40, "r": 20, "t": 40, "b": 40})
        cards.append(
            dbc.Col(
                dbc.Card(
                    dbc.CardBody(
                        [
                            html.Div(
                                opt["label"], className="fw-bold mb-1", style={"fontSize": "12px"}
                            ),
                            dcc.Graph(
                                figure=fig,
                                config={"displayModeBar": False},
                                style={"height": "280px"},
                            ),
                        ]
                    ),
                    className="mb-3",
                ),
                width=6,
            )
        )
    return dbc.Row(cards)
