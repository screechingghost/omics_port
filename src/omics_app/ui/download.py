"""
Port of R tabPanel("💾 Download", ...) (app_12-02.R, lines 1507-1625).

UI-ONLY, as requested -- every downloadButton below is a plain
dbc.Button with no dcc.Download and no callback. None of them fetch
data, build a ZIP, or write to disk. Wiring these up (a callback per
button reading store-dea-results / store-enrichment-results and
returning file bytes, the way ui/visualization.py's TIFF download
buttons work) is future work.

Layout note: R's original is a flat 6-across row of equal-weight
feature boxes (line 1519-1586) plus three more full-width cards
stacked below -- 4 separate blocks competing for attention with no
hierarchy. This restructures the same content, not new content:
  - "Download Everything" is pulled out as one prominent action
    (it's the one most people actually want) instead of being just
    another box in the row.
  - The remaining 5 downloads move to a 3-column grid (was 6 columns
    at width=2, cramped) so each card has room to breathe.
  - Export Format Options collapses from a full card into a single
    inline settings row, since it's a small global toggle, not
    content of its own.
  - R's "Latest Folder Export" status card (line 1588-1598) is wrapped
    in a conditionalPanel that's hidden until a folder export actually
    happens -- since it has nothing to show yet in this UI-only pass
    anyway, it's left out here rather than shown permanently empty;
    add it back (conditionally rendered) once Save to Folder is wired up.
"""

# pyright: reportCallIssue=false, reportInvalidTypeForm=false
import dash_bootstrap_components as dbc
from dash import html


def _download_card(
    icon_title: str, button_label: str, description: str, color: str, accent: str
) -> dbc.Col:
    return dbc.Col(
        dbc.Card(
            dbc.CardBody(
                [
                    html.H6(icon_title, className="mb-3"),
                    html.P(
                        description,
                        style={"fontSize": "12.5px", "color": "var(--muted)", "minHeight": "34px"},
                    ),
                    dbc.Button(
                        button_label, color=color, size="sm", className="w-100 mt-1", disabled=True
                    ),
                ]
            ),
            className=f"h-100 card-accent-{accent}",
        ),
        width=4,
        className="mb-3",
    )


def layout() -> html.Div:
    header = dbc.Card(
        dbc.CardBody(
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H4("Download Results", className="mb-1"),
                            html.P(
                                "Export analysis results, plots, and reports in an organized folder structure.",
                                style={
                                    "fontSize": "13px",
                                    "color": "var(--muted)",
                                    "marginBottom": 0,
                                },
                            ),
                        ],
                        width="auto",
                        className="me-auto",
                    ),
                    dbc.Col(
                        [
                            html.Div(
                                "Plot formats to include",
                                style={
                                    "fontSize": "10.5px",
                                    "fontWeight": "700",
                                    "color": "var(--muted)",
                                    "textTransform": "uppercase",
                                    "letterSpacing": "0.04em",
                                    "marginBottom": "4px",
                                },
                            ),
                            dbc.Checklist(
                                id="export-formats",
                                options=[
                                    {"label": "TIFF (high-res)", "value": "tiff"},
                                    {"label": "SVG (vector)", "value": "svg"},
                                    {"label": "PNG (web)", "value": "png"},
                                ],
                                value=["tiff", "svg"],
                                inline=True,
                                switch=True,
                            ),
                        ],
                        width="auto",
                    ),
                ],
                align="center",
                justify="between",
                className="g-3",
            )
        ),
        className="card-accent-indigo mb-3",
    )

    hero = dbc.Card(
        dbc.CardBody(
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H5(
                                "📦 Complete Package", className="mb-1", style={"color": "white"}
                            ),
                            html.P(
                                "Everything in one archive -- data, plots, and enrichment results.",
                                style={"fontSize": "13px", "color": "#C9CEE8", "marginBottom": 0},
                            ),
                        ],
                        width="auto",
                        className="me-auto",
                    ),
                    dbc.Col(
                        dbc.Button(
                            "Download Everything (ZIP)",
                            id="download-all-zip-btn",
                            color="light",
                            size="md",
                            disabled=True,
                        ),
                        width="auto",
                    ),
                ],
                align="center",
                justify="between",
            )
        ),
        style={"backgroundColor": "var(--indigo)"},
        className="mb-3",
    )

    grid = dbc.Row(
        [
            _download_card(
                "📊 Analysis Data",
                "Download DEA Results",
                "All differential expression analysis results.",
                "success",
                "teal",
            ),
            _download_card(
                "📈 Plots",
                "Download All Plots",
                "All visualizations, in the formats selected above.",
                "info",
                "indigo",
            ),
            _download_card(
                "📄 Report",
                "Download Report",
                "Ready-to-send PDF with separate plot pages.",
                "danger",
                "coral",
            ),
            _download_card(
                "🧬 Enrichment",
                "Download Enrichment Results",
                "GO terms, pathways, and enrichment plots.",
                "warning",
                "amber",
            ),
            _download_card(
                "📁 Save to Folder",
                "Save Everything to Folder",
                "Creates an unzipped export folder on this PC.",
                "secondary",
                "indigo",
            ),
        ],
        className="g-3",
    )

    individual_downloads = dbc.Card(
        dbc.CardBody(
            [
                html.H5("Individual Downloads", className="mb-1"),
                html.P(
                    "Download a specific comparison or analysis on its own.",
                    style={"fontSize": "13px", "color": "var(--muted)"},
                ),
                html.Div(
                    "Run an analysis to see per-comparison downloads here.",
                    style={
                        "color": "var(--muted)",
                        "fontSize": "13px",
                        "padding": "20px",
                        "textAlign": "center",
                        "border": "1px dashed var(--border)",
                        "borderRadius": "8px",
                    },
                ),
            ]
        ),
        className="mt-1",
    )

    return html.Div([header, hero, grid, individual_downloads])
