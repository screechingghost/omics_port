"""
Port of R tabPanel("📁 Data Upload", ...) (app_12-02.R, lines 904-963).

Server-side logic this replaces: the file-load reactive chain that
reads the uploaded Excel file, lets the user pick a sheet and row-name
column, then calls build_main_data_index() to auto-classify columns.

Performance note: uploaded file bytes are cached server-side (see
omics_app.server.cache) and referenced by a small token in dcc.Store,
rather than shuttling the full base64-encoded file through the browser
on every callback. The earlier version stored raw base64 content
directly in a client-side dcc.Store, which meant the whole file got
re-serialized/re-sent on every step (upload -> sheet select -> load
click) -- the actual cause of slow-feeling uploads, worse the larger
the file.
"""

import base64
import io
import uuid

import dash_bootstrap_components as dbc
import pandas as pd
from dash import Input, Output, State, callback, ctx, dash_table, dcc, html

from omics_app.data.columns import build_main_data_index
from omics_app.server import cache


def layout() -> html.Div:
    return html.Div(
        [
            html.Div(id="upload-gate-banner", className="mb-3"),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Div(
                                                        "01",
                                                        style={
                                                            "width": "34px",
                                                            "height": "34px",
                                                            "borderRadius": "50%",
                                                            "backgroundColor": "var(--indigo)",
                                                            "color": "white",
                                                            "display": "flex",
                                                            "alignItems": "center",
                                                            "justifyContent": "center",
                                                            "fontWeight": "700",
                                                            "fontSize": "13px",
                                                            "marginRight": "12px",
                                                        },
                                                    ),
                                                    html.Div(
                                                        [
                                                            html.H4(
                                                                "Import your dataset",
                                                                className="mb-1",
                                                                style={"fontSize": "19px"},
                                                            ),
                                                            html.Div(
                                                                "Start by choosing the Excel workbook containing your measurements.",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                ],
                                                style={
                                                    "display": "flex",
                                                    "alignItems": "center",
                                                    "marginBottom": "18px",
                                                },
                                            ),
                                            html.Div(
                                                [
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                "File type",
                                                                style={
                                                                    "fontSize": "11px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                            html.Div(
                                                                "Excel workbook",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "fontWeight": "600",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                "Accepted formats",
                                                                style={
                                                                    "fontSize": "11px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                            html.Div(
                                                                ".xlsx  ·  .xls",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "fontWeight": "600",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                "Contains",
                                                                style={
                                                                    "fontSize": "11px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                            html.Div(
                                                                "Measurements + metadata",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "fontWeight": "600",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                ],
                                                style={
                                                    "display": "flex",
                                                    "justifyContent": "space-between",
                                                    "gap": "12px",
                                                    "padding": "12px 14px",
                                                    "backgroundColor": "#F8FAFC",
                                                    "border": "1px solid var(--border)",
                                                    "borderRadius": "8px",
                                                    "marginBottom": "16px",
                                                },
                                            ),
                                            dcc.Upload(
                                                id="main-file-upload",
                                                children=html.Div(
                                                    [
                                                        html.Div(
                                                            "↑",
                                                            style={
                                                                "width": "42px",
                                                                "height": "42px",
                                                                "borderRadius": "50%",
                                                                "backgroundColor": "#EDE9FE",
                                                                "color": "var(--indigo)",
                                                                "display": "flex",
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
                                                                "fontSize": "25px",
                                                                "fontWeight": "500",
                                                                "margin": "0 auto 10px",
                                                            },
                                                        ),
                                                        html.Div(
                                                            [
                                                                html.Strong(
                                                                    "Choose an Excel workbook"
                                                                ),
                                                                html.Span(
                                                                    " or drag it here",
                                                                    style={"color": "var(--muted)"},
                                                                ),
                                                            ],
                                                            style={"fontSize": "14px"},
                                                        ),
                                                        html.Div(
                                                            "Your file remains available while you configure sheets and columns.",
                                                            style={
                                                                "fontSize": "12px",
                                                                "color": "var(--muted)",
                                                                "marginTop": "5px",
                                                            },
                                                        ),
                                                    ]
                                                ),
                                                accept=".xlsx,.xls",
                                                style={
                                                    "border": "1.5px dashed var(--indigo)",
                                                    "borderRadius": "10px",
                                                    "padding": "22px 16px",
                                                    "textAlign": "center",
                                                    "backgroundColor": "#FAFAFF",
                                                    "cursor": "pointer",
                                                },
                                            ),
                                            html.Div(id="main-file-name-display", className="mt-3"),
                                            html.Div(id="main-upload-status", className="mt-2"),
                                            dcc.Upload(
                                                id="replace-file-upload",
                                                children=dbc.Button(
                                                    ["↻ ", "Choose a different workbook"],
                                                    color="outline-primary",
                                                    size="sm",
                                                ),
                                                accept=".xlsx,.xls",
                                                style={"display": "none", "marginTop": "12px"},
                                            ),
                                        ]
                                    ),
                                    html.Div(
                                        [
                                            dbc.Label("Select Sheet:"),
                                            dcc.Dropdown(
                                                id="sheet-selector",
                                                options=[],
                                                value=None,
                                                clearable=False,
                                            ),
                                        ],
                                        className="mt-3",
                                    ),
                                    html.Div(
                                        [
                                            dbc.Label("Row Name Column:"),
                                            dcc.Dropdown(
                                                id="rowname-selector",
                                                options=[],
                                                value=None,
                                                clearable=False,
                                            ),
                                        ],
                                        className="mt-3",
                                    ),
                                    dbc.Button(
                                        "Load Data",
                                        id="load-data-btn",
                                        color="primary",
                                        className="mt-3 w-100",
                                    ),
                                ]
                            ),
                            className="card-accent-indigo",
                        ),
                        width=6,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Div(
                                                        "02",
                                                        style={
                                                            "width": "34px",
                                                            "height": "34px",
                                                            "borderRadius": "50%",
                                                            "backgroundColor": "var(--teal)",
                                                            "color": "white",
                                                            "display": "flex",
                                                            "alignItems": "center",
                                                            "justifyContent": "center",
                                                            "fontWeight": "700",
                                                            "fontSize": "13px",
                                                            "marginRight": "12px",
                                                        },
                                                    ),
                                                    html.Div(
                                                        [
                                                            html.H4(
                                                                "Import existing results",
                                                                className="mb-1",
                                                                style={"fontSize": "19px"},
                                                            ),
                                                            html.Div(
                                                                "Optional — use a previous analysis for visualization without rerunning it.",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                ],
                                                style={
                                                    "display": "flex",
                                                    "alignItems": "center",
                                                    "marginBottom": "18px",
                                                },
                                            ),
                                            html.Div(
                                                [
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                "File type",
                                                                style={
                                                                    "fontSize": "11px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                            html.Div(
                                                                "Delimited results file",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "fontWeight": "600",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                "Accepted formats",
                                                                style={
                                                                    "fontSize": "11px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                            html.Div(
                                                                ".txt  ·  .csv",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "fontWeight": "600",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                "Used for",
                                                                style={
                                                                    "fontSize": "11px",
                                                                    "color": "var(--muted)",
                                                                },
                                                            ),
                                                            html.Div(
                                                                "Plots + exploration",
                                                                style={
                                                                    "fontSize": "13px",
                                                                    "fontWeight": "600",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                ],
                                                style={
                                                    "display": "flex",
                                                    "justifyContent": "space-between",
                                                    "gap": "12px",
                                                    "padding": "12px 14px",
                                                    "backgroundColor": "#F4FBFA",
                                                    "border": "1px solid var(--border)",
                                                    "borderRadius": "8px",
                                                    "marginBottom": "16px",
                                                },
                                            ),
                                            dcc.Upload(
                                                id="sig-results-upload",
                                                children=html.Div(
                                                    [
                                                        html.Div(
                                                            "↑",
                                                            style={
                                                                "width": "42px",
                                                                "height": "42px",
                                                                "borderRadius": "50%",
                                                                "backgroundColor": "#DDF7F4",
                                                                "color": "var(--teal)",
                                                                "display": "flex",
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
                                                                "fontSize": "25px",
                                                                "fontWeight": "500",
                                                                "margin": "0 auto 10px",
                                                            },
                                                        ),
                                                        html.Div(
                                                            [
                                                                html.Strong(
                                                                    "Choose a results file"
                                                                ),
                                                                html.Span(
                                                                    " or drag it here",
                                                                    style={"color": "var(--muted)"},
                                                                ),
                                                            ],
                                                            style={"fontSize": "14px"},
                                                        ),
                                                        html.Div(
                                                            "Select the result type below, then load it for visualization.",
                                                            style={
                                                                "fontSize": "12px",
                                                                "color": "var(--muted)",
                                                                "marginTop": "5px",
                                                            },
                                                        ),
                                                    ]
                                                ),
                                                accept=".txt,.csv",
                                                style={
                                                    "border": "1.5px dashed var(--teal)",
                                                    "borderRadius": "10px",
                                                    "padding": "22px 16px",
                                                    "textAlign": "center",
                                                    "backgroundColor": "#FAFFFE",
                                                    "cursor": "pointer",
                                                },
                                            ),
                                            html.Div(
                                                [
                                                    dbc.Label(
                                                        "Results file contains:",
                                                        style={
                                                            "fontSize": "13px",
                                                            "fontWeight": "600",
                                                        },
                                                    ),
                                                    dcc.Dropdown(
                                                        id="sig-results-type",
                                                        options=[
                                                            {
                                                                "label": "DEA Results",
                                                                "value": "dea",
                                                            },
                                                            {
                                                                "label": "Significant Proteins (p-value)",
                                                                "value": "sig_p",
                                                            },
                                                            {
                                                                "label": "Significant Proteins (q-value)",
                                                                "value": "sig_q",
                                                            },
                                                        ],
                                                        value="dea",
                                                        clearable=False,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Button(
                                                "Load Significant Results",
                                                id="load-sig-btn",
                                                color="success",
                                                className="mt-3 w-100",
                                            ),
                                            html.Div(id="sig-load-status", className="mt-2"),
                                        ]
                                    ),
                                ]
                            ),
                            className="card-accent-teal",
                        ),
                        width=6,
                    ),
                ],
                className="mb-3",
            ),
            dbc.Row(
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4("Data Preview"),
                                dcc.Loading(html.Div(id="data-preview-table")),
                            ]
                        )
                    ),
                    width=12,
                ),
                className="mb-3",
            ),
            dbc.Row(
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4("Dataset Information"),
                                html.Div(id="data-info-text"),
                            ]
                        )
                    ),
                    width=12,
                )
            ),
            # Small token only -- the actual decoded file bytes live in the
            # server-side cache (see omics_app.server.cache), not here.
            dcc.Store(id="store-uploaded-file-content"),
        ]
    )


def _decode_upload(contents: str) -> bytes:
    _content_type, content_string = contents.split(",", 1)
    return base64.b64decode(content_string)


def _format_file_size(size_bytes: int) -> str:
    if size_bytes < 1024:
        return f"{size_bytes} B"
    if size_bytes < 1024**2:
        return f"{size_bytes / 1024:.1f} KB"
    return f"{size_bytes / 1024**2:.2f} MB"


def _read_excel_fast(file_like, **kwargs):
    """
    python-calamine reads Excel files roughly 20x faster than the
    default openpyxl engine (measured: 278ms vs 14ms on a 500-row test
    file, and the gap widens on larger files) since it's a Rust parser
    rather than pure Python. Falls back to openpyxl automatically if
    python-calamine isn't installed, so this never hard-fails.
    """
    try:
        return pd.read_excel(file_like, engine="calamine", **kwargs)
    except ImportError:
        file_like.seek(0)
        return pd.read_excel(file_like, engine="openpyxl", **kwargs)


# --- Gate: require user name + ID (set on the Home tab) before allowing
# any upload action, mirroring the R app's "User information required"
# warning on Home, but enforced here where it actually blocks progress.


def _has_user_info(app_config: dict | None) -> bool:
    app_config = app_config or {}
    return bool((app_config.get("user_name") or "").strip()) and bool(
        (app_config.get("user_id") or "").strip()
    )


@callback(
    Output("upload-gate-banner", "children"),
    Output("load-data-btn", "disabled"),
    Output("load-sig-btn", "disabled"),
    Input("store-app-config", "data"),
)
def gate_upload_on_user_info(app_config):
    if _has_user_info(app_config):
        return None, False, False
    banner = dbc.Alert(
        [
            html.Strong("User Name and User ID required. "),
            "Please fill them in on the Home tab before loading data -- "
            "they're used to organize exported results into user-specific folders.",
        ],
        color="warning",
    )
    return banner, True, True


@callback(
    Output("store-uploaded-file-content", "data"),
    Output("sheet-selector", "options"),
    Output("sheet-selector", "value"),
    Output("main-file-name-display", "children"),
    Output("main-upload-status", "children"),
    Output("main-file-upload", "style"),
    Output("replace-file-upload", "style"),
    Input("main-file-upload", "contents"),
    Input("replace-file-upload", "contents"),
    State("main-file-upload", "filename"),
    State("replace-file-upload", "filename"),
)
def on_file_uploaded(main_contents, replace_contents, main_filename, replace_filename):
    """
    Reads sheet names as soon as a file is dropped. The decoded bytes go
    into the server-side cache under a short-lived token; only that
    token (plus the filename) goes into the client-side dcc.Store, so
    the large file content never round-trips through the browser again
    on subsequent callbacks (sheet select, load click).
    """
    triggered_id = ctx.triggered_id

    if triggered_id == "replace-file-upload":
        contents = replace_contents
        filename = replace_filename
    else:
        contents = main_contents
        filename = main_filename

    if contents is None:
        return (
            None,
            [],
            None,
            None,
            None,
            {
                "display": "block",
                "border": "1.5px dashed var(--indigo)",
                "borderRadius": "10px",
                "padding": "22px 16px",
                "textAlign": "center",
                "backgroundColor": "#FAFAFF",
                "cursor": "pointer",
            },
            {"display": "none"},
        )
    decoded = _decode_upload(contents)
    try:
        sheet_names = pd.ExcelFile(io.BytesIO(decoded), engine="calamine").sheet_names
    except ImportError:
        sheet_names = pd.ExcelFile(io.BytesIO(decoded), engine="openpyxl").sheet_names

    token = str(uuid.uuid4())
    cache.set(f"upload:{token}", decoded)

    options = [{"label": s, "value": s} for s in sheet_names]
    default_value = sheet_names[0] if sheet_names else None
    file_size = _format_file_size(len(decoded))

    filename_display = dbc.Card(
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.Span("📄", style={"fontSize": "22px", "marginRight": "10px"}),
                        html.Div(
                            [
                                html.Div(
                                    filename or "Uploaded workbook",
                                    style={
                                        "fontWeight": "600",
                                        "color": "var(--ink)",
                                        "wordBreak": "break-word",
                                    },
                                ),
                                html.Small(
                                    f"Excel workbook · {file_size}",
                                    style={"color": "var(--muted)"},
                                ),
                            ]
                        ),
                    ],
                    style={"display": "flex", "alignItems": "center"},
                ),
            ]
        ),
        className="border-0 shadow-sm",
        style={"backgroundColor": "#F8FAFC"},
    )

    upload_status = dbc.Alert(
        [
            html.Strong("Upload successful"),
            html.Span(f" · {len(sheet_names)} sheet(s) detected and ready to configure."),
        ],
        color="success",
        className="py-2 mb-0",
    )

    return (
        {"filename": filename, "token": token},
        options,
        default_value,
        filename_display,
        upload_status,
        {"display": "none"},
        {"display": "inline-block", "marginTop": "12px"},
    )


@callback(
    Output("rowname-selector", "options"),
    Output("rowname-selector", "value"),
    Input("sheet-selector", "value"),
    State("store-uploaded-file-content", "data"),
    prevent_initial_call=True,
)
def on_sheet_selected(sheet_name, stored_file):
    """Reads column headers for the chosen sheet and populates
    rowname-selector's options (component always exists -- see note in
    on_file_uploaded above)."""
    if not stored_file or not sheet_name:
        return [], None

    decoded = cache.get(f"upload:{stored_file['token']}")
    if decoded is None:
        return [], None  # cache expired -- user should re-upload

    header_only = _read_excel_fast(io.BytesIO(decoded), sheet_name=sheet_name, nrows=0)
    columns = list(header_only.columns)

    options = [{"label": c, "value": c} for c in columns]
    return options, (columns[0] if columns else None)


def _column_badge_row(label: str, columns: list[str], accent: str) -> html.Div:
    """One labeled row of badge chips -- used for group/abundance/metadata
    column lists in the Dataset Information panel. Chips scan much
    faster than a long comma-separated sentence, especially once a
    real file has 20+ abundance columns."""
    if not columns:
        chips = html.Span("none detected", style={"color": "var(--muted)", "fontSize": "13px"})
    else:
        chips = html.Div(
            [
                dbc.Badge(
                    col,
                    color=None,
                    className="me-1 mb-1",
                    style={
                        "backgroundColor": f"var(--{accent})",
                        "color": "#FFFFFF",
                        "fontWeight": "500",
                        "fontSize": "12px",
                    },
                )
                for col in columns
            ],
            style={"display": "flex", "flexWrap": "wrap", "gap": "2px"},
        )
    return html.Div(
        [
            html.Div(
                f"{label} ({len(columns)})",
                style={"fontWeight": "600", "fontSize": "13px", "marginBottom": "6px"},
            ),
            chips,
        ],
        className="mb-3",
    )


def _info_stat(label: str, value, accent: str) -> dbc.Col:
    """One stat block in the Dataset Information panel -- colored to
    match the accent system used elsewhere (indigo/teal/coral/amber)."""
    return dbc.Col(
        html.Div(
            [
                html.Div(
                    str(value),
                    style={
                        "fontSize": "28px",
                        "fontWeight": "700",
                        "color": f"var(--{accent})",
                    },
                ),
                html.Div(
                    label,
                    style={
                        "fontSize": "13px",
                        "color": "var(--muted)",
                        "marginTop": "2px",
                    },
                ),
            ],
            style={"textAlign": "center", "padding": "12px 0"},
        ),
        width=3,
    )


@callback(
    Output("store-main-data", "data"),
    Output("data-preview-table", "children"),
    Output("data-info-text", "children"),
    Output("store-comparisons", "data", allow_duplicate=True),
    Input("load-data-btn", "n_clicks"),
    State("store-uploaded-file-content", "data"),
    State("sheet-selector", "value"),
    State("rowname-selector", "value"),
    prevent_initial_call=True,
)
def on_load_data_clicked(n_clicks, stored_file, sheet_name, rowname_col):
    """
    Port of the "Load Data" action -- reads the full sheet, runs
    build_main_data_index() for column auto-detection, and populates
    the preview table + dataset info panel.

    Also clears store-comparisons (Output #4) on every successful load.
    Any comparisons configured against a *previous* file reference that
    file's group/abundance column names -- silently leaving them in
    place after a new file loads is what let stale dropdown values
    reach Comparisons' rebuilt dropdowns with a new, non-matching
    options list (see the note in comparisons.py's _comparison_card
    for the resulting "Cannot read properties of null" crash this
    caused in the browser).
    """
    if not stored_file or not sheet_name:
        return None, None, dbc.Alert("No file loaded yet.", color="secondary"), None

    decoded = cache.get(f"upload:{stored_file['token']}")
    if decoded is None:
        return (
            None,
            None,
            dbc.Alert("Upload session expired -- please re-upload the file.", color="warning"),
            None,
        )

    df = _read_excel_fast(io.BytesIO(decoded), sheet_name=sheet_name)

    if rowname_col and rowname_col in df.columns:
        df = df.set_index(rowname_col)

    display_df = df.reset_index()
    column_index = build_main_data_index(display_df)

    preview = dash_table.DataTable(
        data=display_df.head(100).to_dict(
            "records"
        ),  # matches R's head(rv$main_data, 100), line 2014
        columns=[{"name": c, "id": c} for c in display_df.columns],
        page_size=10,
        style_table={"overflowX": "auto"},
        style_as_list_view=False,
        style_cell={
            "fontSize": "12px",
            "fontFamily": "monospace",
            "textAlign": "left",
            "padding": "8px",
            "border": "1px solid #D8DCE8",
        },
        style_header={
            "textAlign": "center",
            "fontWeight": "600",
            "backgroundColor": "var(--surface-alt, #F7F8FC)",
            "border": "1px solid #D8DCE8",
        },
    )

    info_panel = html.Div(
        [
            dbc.Row(
                [
                    _info_stat("Rows", len(df), "indigo"),
                    _info_stat("Columns", len(display_df.columns), "teal"),
                    _info_stat("Group Columns", len(column_index["group_cols"]), "coral"),
                    _info_stat("Abundance Columns", len(column_index["abundance_cols"]), "amber"),
                ],
            ),
            html.Hr(),
            dbc.Row(
                [
                    dbc.Col(
                        html.P([html.Strong("File: "), stored_file["filename"]], className="mb-1"),
                        width=6,
                    ),
                    dbc.Col(
                        html.P([html.Strong("Sheet: "), sheet_name], className="mb-1"),
                        width=6,
                    ),
                ],
                className="mb-2",
            ),
            _column_badge_row("Group Columns", column_index["group_cols"], "coral"),
            _column_badge_row("Abundance Columns", column_index["abundance_cols"], "amber"),
            _column_badge_row("Metadata Candidates", column_index["metadata_candidates"], "indigo"),
        ]
    )

    # Stored for downstream tabs (Comparisons, Analysis, etc.) to consume.
    store_payload = {
        "filename": stored_file["filename"],
        "sheet_name": sheet_name,
        "rowname_col": rowname_col,
        "data": display_df.to_dict("records"),
        "column_index": column_index,
    }
    return store_payload, preview, info_panel, None


@callback(
    Output("store-dea-results", "data", allow_duplicate=True),
    Output("sig-load-status", "children"),
    Input("load-sig-btn", "n_clicks"),
    State("sig-results-upload", "contents"),
    State("sig-results-upload", "filename"),
    State("store-dea-results", "data"),
    prevent_initial_call=True,
)
def on_load_sig_results_clicked(n_clicks, contents, filename, existing_dea_results):
    """
    Port of observeEvent(input$load_sig_btn, ...) (R line 1986). Loads a
    pre-existing tab-separated significant-results file for
    visualization-only use, without running a new DEA.
    """
    if contents is None:
        return existing_dea_results, dbc.Alert(
            "Please choose a file before clicking Load.", color="warning"
        )

    try:
        decoded = _decode_upload(contents)
        sig_data = pd.read_csv(io.BytesIO(decoded), sep="\t")
    except Exception as e:
        return existing_dea_results, dbc.Alert(
            f"Error loading significant results: {e}", color="danger"
        )

    dea_results = dict(existing_dea_results or {})
    dea_results["uploaded_sig"] = {
        "significant_pvalue": sig_data.to_dict("records"),
        "all_proteins": sig_data.to_dict("records"),
        "source_filename": filename,
    }
    status = dbc.Alert(
        f"Significant results loaded successfully! ({len(sig_data)} proteins from {filename})",
        color="success",
    )
    return dea_results, status
