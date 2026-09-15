"""
Port of R tabPanel("ℹ️ Session Info", ...) (app_12-02.R, lines 1630-1659)
and its three renderers:
  output$session_info_text  (line 15241) -- live snapshot of app state
  output$analysis_log_text  (line 15256) -- rv$analysis_log, joined
  output$packages_info_text (line 15260) -- R version + loaded packages

What's wired:
  - Session Information: a live callback reading every dcc.Store in the
    app (store-app-config, store-main-data, store-comparisons,
    store-dea-results, store-enrichment-results, store-color-mapping).
  - Analysis Log: renders store-analysis-log, which upload.py,
    comparisons.py, analysis.py, ui/enrichment.py, and colors.py each
    append a line to on their respective successful actions (see
    stats/session_log.py). Each entry carries a timestamp -- a
    deliberate small improvement over R's plain message list, not a
    literal port requirement.
  - Loaded Packages: real installed package versions via
    importlib.metadata, computed once when this module's layout() runs
    at app startup.
"""

# pyright: reportCallIssue=false, reportInvalidTypeForm=false
import sys
from importlib import metadata as importlib_metadata

import dash_bootstrap_components as dbc
from dash import Input, Output, callback, html

_METHOD_LABELS = {"normal": "T-Test", "anova": "ANOVA (Multi-Group)"}

_PACKAGES_OF_INTEREST = [
    "dash",
    "dash-bootstrap-components",
    "dash-daq",
    "flask-caching",
    "pandas",
    "numpy",
    "scipy",
    "scikit-learn",
    "statsmodels",
    "plotly",
    "kaleido",
    "Pillow",
    "gseapy",
]


def _package_versions() -> list[tuple[str, str]]:
    versions = []
    for name in _PACKAGES_OF_INTEREST:
        try:
            versions.append((name, importlib_metadata.version(name)))
        except importlib_metadata.PackageNotFoundError:
            versions.append((name, "not installed"))
    return versions


def _info_card(title: str, body_id: str, placeholder: str) -> dbc.Card:
    return dbc.Card(
        dbc.CardBody(
            [
                html.H3(title, style={"color": "var(--ink)"}),
                html.Pre(
                    placeholder,
                    id=body_id,
                    style={
                        "fontSize": "12.5px",
                        "color": "var(--muted)",
                        "backgroundColor": "var(--bg)",
                        "padding": "12px",
                        "borderRadius": "8px",
                        "whiteSpace": "pre-wrap",
                        "marginBottom": 0,
                    },
                ),
            ]
        ),
        className="mb-3",
    )


def layout() -> html.Div:
    packages_lines = [f"Python Version: {sys.version.split()[0]}", "", "Key Packages:"]
    packages_lines += [f"  {name}: {version}" for name, version in _package_versions()]

    return html.Div(
        [
            _info_card("Session Information", "session-info-text", "Loading session info..."),
            _info_card("Analysis Log", "analysis-log-text", "No actions logged yet."),
            dbc.Card(
                dbc.CardBody(
                    [
                        html.H3("Loaded Packages", style={"color": "var(--ink)"}),
                        html.Pre(
                            "\n".join(packages_lines),
                            style={
                                "fontSize": "12.5px",
                                "color": "var(--muted)",
                                "backgroundColor": "var(--bg)",
                                "padding": "12px",
                                "borderRadius": "8px",
                                "whiteSpace": "pre-wrap",
                                "marginBottom": 0,
                            },
                        ),
                    ]
                ),
                className="mb-3",
            ),
        ]
    )


@callback(
    Output("session-info-text", "children"),
    Input("store-app-config", "data"),
    Input("store-main-data", "data"),
    Input("store-comparisons", "data"),
    Input("store-dea-results", "data"),
    Input("store-enrichment-results", "data"),
    Input("store-color-mapping", "data"),
)
def render_session_info(
    app_config, main_data, comparisons_data, dea_results, enrichment_results, color_mapping
):
    """Port of output$session_info_text (R line 15241)."""
    app_config = app_config or {}
    analysis_type = app_config.get("analysis_type")
    comparison_method = app_config.get("comparison_method")
    user_name = (app_config.get("user_name") or "").strip()
    user_id = (app_config.get("user_id") or "").strip()

    lines = [
        f"Analysis Type: {analysis_type.title() if analysis_type else 'Not set'}",
        f"Comparison Method: {_METHOD_LABELS.get(comparison_method, comparison_method or 'Not set')}",
        f"User Name: {user_name or 'Not set'}",
        f"User ID: {user_id or 'Not set'}",
        f"Input File: {(main_data or {}).get('filename', 'Not recorded')}",
        f"Input Sheet: {(main_data or {}).get('sheet_name', 'Not recorded')}",
        f"Row Name Column: {(main_data or {}).get('rowname_col', 'Not recorded')}",
        f"Data Loaded: {main_data is not None}",
        f"Comparisons Created: {len((comparisons_data or {}).get('comparisons') or [])}",
        f"Analyses Completed: {len(dea_results or {})}",
        f"Enrichment Results: {len(enrichment_results or {})}",
        f"Colors Mapped: {len(color_mapping or {})}",
    ]
    return "\n".join(lines)


@callback(
    Output("analysis-log-text", "children"),
    Input("store-analysis-log", "data"),
)
def render_analysis_log(log_entries):
    """Port of output$analysis_log_text (R line 15256)."""
    if not log_entries:
        return "No actions logged yet."
    return "\n".join(log_entries)
