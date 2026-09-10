"""
Port of R tabPanel("🧬 Enrichment", ...) (app_12-02.R, lines 1385-1502).

UI-ONLY, as requested -- this module is layout, no callbacks. Nothing
here reads store-dea-results, calls any enrichment library, or writes
to a store; every control keeps R's default value but is otherwise
inert. Wiring it up (comparison selector driven by store-dea-results,
Run Enrichment Analysis -> a real GO/Reactome or DAVID call, results
table/plot rendering) is future work, not part of this pass -- follow
ui/analysis.py's pattern (@callback functions reading/writing
dcc.Store) when that's ready.

Known simplification even at the UI level: R's two engine-specific
settings blocks (Bioconductor GO/Reactome vs. DAVID Webservice, R
lines 1404-1439) are shown via conditionalPanel, toggled by the
enrichment_engine radio. Reproducing that here would need a callback
(even a purely client-side visibility toggle is still a callback), so
for this UI-only pass both blocks are simply shown at once, each
clearly labeled. Add a toggle callback for `enrichment-engine` when
this tab gets wired up.

Organism list (line 272-283 in R, get_organism_options()) is
abbreviated to a representative subset here -- the full mapping needs
that function ported too, which is backend work, not layout.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html

ORGANISM_OPTIONS = [
    {"label": "Human (GO + Reactome)", "value": "human"},
    {"label": "Mouse (GO + Reactome)", "value": "mouse"},
    {"label": "Rat (GO + Reactome)", "value": "rat"},
    {"label": "Zebrafish (GO + Reactome)", "value": "zebrafish"},
    {"label": "Fly (GO only)", "value": "fly"},
    {"label": "Worm (C. elegans) (GO only)", "value": "worm"},
    {"label": "Yeast (GO + Reactome)", "value": "yeast"},
    {"label": "Arabidopsis (GO only)", "value": "arabidopsis"},
]

PLOT_TYPE_OPTIONS = [
    {"label": "GO: By Ontology (BP/MF/CC)", "value": "go_ontology"},
    {"label": "Grouped Barplot", "value": "barplot"},
    {"label": "Dotplot", "value": "dotplot"},
    {"label": "Bubble Plot", "value": "bubble"},
    {"label": "Up vs Down Comparison", "value": "up_down_compare"},
    {"label": "Gene-Pathway Network (cnetplot)", "value": "cnetplot"},
    {"label": "Enrichment Map (term similarity)", "value": "emapplot"},
    {"label": "Chord: Pathway Up", "value": "chord_up"},
    {"label": "Chord: Pathway Down", "value": "chord_down"},
]


def layout() -> html.Div:
    settings = dbc.Card(
        dbc.CardBody(
            [
                html.H4("Enrichment Analysis Setup"),
                dbc.Label("Comparison:", className="mt-2"),
                dcc.Dropdown(id="enrichment-comparison-select", options=[], placeholder="Run Analysis first"),

                dbc.Label("Select Organism:", className="mt-3"),
                dcc.Dropdown(id="enrichment-organism", options=ORGANISM_OPTIONS, value="human", clearable=False),

                dbc.Label("Enrichment Engine:", className="mt-3"),
                dbc.RadioItems(
                    id="enrichment-engine",
                    options=[
                        {"label": "Bioconductor GO/Reactome", "value": "bioc"},
                        {"label": "DAVID Webservice", "value": "david"},
                    ],
                    value="bioc",
                ),

                # Both engine-specific blocks are always shown -- see module
                # docstring (no conditionalPanel-style toggle without a callback).
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
                        html.H6("DAVID settings", style={"color": "var(--muted)"}),
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
                            value="OFFICIAL_GENE_SYMBOL", clearable=False,
                        ),
                        dbc.Label("DAVID Categories:", className="mt-2"),
                        dbc.Input(id="david-categories", type="text",
                                   value="GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,KEGG_PATHWAY"),
                        dbc.Label("DAVID EASE threshold:", className="mt-2"),
                        dbc.Input(id="david-ease-threshold", type="number", value=1, min=0.001, max=1, step=0.005),
                        dbc.Label("DAVID minimum count:", className="mt-2"),
                        dbc.Input(id="david-count-threshold", type="number", value=1, min=1, max=100, step=1),
                        dbc.Checklist(
                            id="david-mock",
                            options=[{"label": "Use mock DAVID for local testing only", "value": "mock"}],
                            value=[], switch=True, className="mt-2",
                        ),
                    ]
                ),

                html.Hr(),
                dbc.Label("Enrichment p-value cutoff:"),
                dbc.Input(id="enrichment-pvalue-cutoff", type="number", value=0.05, min=0.001, max=1, step=0.005),
                dbc.Label("Enrichment q-value cutoff:", className="mt-2"),
                dbc.Input(id="enrichment-qvalue-cutoff", type="number", value=1, min=0.001, max=1, step=0.005),

                dbc.Button("🧬 Run Enrichment Analysis", id="run-enrichment-btn",
                           color="primary", size="lg", className="w-100 mt-3", disabled=True),
                dbc.Button("▶ Run All Enrichment Analyses", id="run-all-enrichment-btn",
                           color="info", className="w-100 mt-2", disabled=True),

                html.Hr(),
                html.H5("Progress:"),
                html.Pre("Not run yet.", id="enrichment-progress",
                          style={"fontSize": "12px", "color": "var(--muted)", "backgroundColor": "var(--bg)",
                                 "padding": "8px", "borderRadius": "6px"}),
            ]
        ),
        className="card-accent-teal",
    )

    results = dbc.Card(
        dbc.CardBody(
            [
                html.H4("Enrichment Results"),
                html.Pre("Run an enrichment analysis to see a summary here.", id="enrichment-summary",
                          style={"fontSize": "13px", "color": "var(--muted)"}),
                html.Hr(),
                html.H5("Enrichment Table:"),
                html.Div(
                    "No results yet.",
                    id="enrichment-results-table",
                    style={"color": "var(--muted)", "fontSize": "13px", "padding": "24px",
                           "textAlign": "center", "border": "1px dashed var(--border)", "borderRadius": "8px"},
                ),
                html.Hr(),
                dbc.Row(
                    [
                        dbc.Col(html.H5("Enrichment Visualization:"), width=6),
                        dbc.Col(
                            dcc.Dropdown(id="enrichment-plot-type", options=PLOT_TYPE_OPTIONS,
                                         value="go_ontology", clearable=False),
                            width=6,
                        ),
                    ]
                ),
                html.Div(
                    "Plot will appear here after running an enrichment analysis.",
                    id="enrichment-plot-placeholder",
                    style={"height": "400px", "display": "flex", "alignItems": "center", "justifyContent": "center",
                           "color": "var(--muted)", "border": "1px dashed var(--border)", "borderRadius": "8px",
                           "marginTop": "12px"},
                ),
                dbc.Button("⬇ Download TIFF", id="download-enrichment-plot-btn",
                           color="success", size="sm", outline=True, className="mt-2", disabled=True),
                html.Hr(),
                html.H5("All Plot Preview"),
                html.Div(id="enrichment-all-plots-preview", style={"color": "var(--muted)", "fontSize": "13px"}),
            ]
        ),
        className="card-accent-indigo",
    )

    return dbc.Row([dbc.Col(settings, width=4), dbc.Col(results, width=8)])
