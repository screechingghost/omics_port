import dash_bootstrap_components as dbc
from dash import Input, Output, callback, dcc, html

ANALYSIS_TYPE_OPTIONS = [
    {"label": "Proteomics", "value": "proteomics"},
    {"label": "Metabolomics", "value": "metabolomics"},
    {"label": "Phosphomatics", "value": "phosphomatics"},
]

COMPARISON_METHOD_OPTIONS = [
    {"label": "T-Test", "value": "normal"},
    {"label": "ANOVA (Multi-Group)", "value": "anova"},
]

_METHOD_LABELS = {"normal": "T-Test", "anova": "ANOVA (Multi-Group)"}


def layout() -> html.Div:
    return html.Div(
        [
            dbc.Row(
                dbc.Col(
                    html.Div(
                        [
                            html.H1(
                                "Welcome to the Multi-Omics Analysis Suite",
                                style={"color": "var(--ink)", "marginBottom": "6px"},
                            ),
                            html.P(
                                "Comprehensive platform for Proteomics, Metabolomics, "
                                "and Phosphomatics analysis",
                                style={"fontSize": "16px", "color": "var(--muted)"},
                            ),
                            html.Div(
                                style={
                                    "height": "3px",
                                    "width": "56px",
                                    "background": "var(--indigo)",
                                    "borderRadius": "2px",
                                    "margin": "16px 0",
                                }
                            ),
                            html.P(
                                "Select your analysis type and comparison method to get started",
                                style={"fontSize": "15px", "color": "var(--ink)"},
                            ),
                        ]
                    ),
                    width=12,
                ),
                className="mb-3",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H4("📊 Analysis Type", style={"color": "var(--ink)"}),
                                    dcc.Dropdown(
                                        id="analysis-type",
                                        options=ANALYSIS_TYPE_OPTIONS,
                                        value="proteomics",
                                        clearable=False,
                                    ),
                                    html.P(
                                        "Choose the type of omics data you want to analyze",
                                        style={"fontSize": "13px", "color": "var(--muted)"},
                                        className="mt-2",
                                    ),
                                ]
                            ),
                            className="card-accent-indigo",
                        ),
                        width=4,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H4("🔬 Comparison Method", style={"color": "var(--ink)"}),
                                    dcc.Dropdown(
                                        id="comparison-method",
                                        options=COMPARISON_METHOD_OPTIONS,
                                        value="normal",
                                        clearable=False,
                                    ),
                                    html.P(
                                        "Select comparison strategy for your experimental design",
                                        style={"fontSize": "13px", "color": "var(--muted)"},
                                        className="mt-2",
                                    ),
                                ]
                            ),
                            className="card-accent-teal",
                        ),
                        width=4,
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H4("💾 Session Management", style={"color": "var(--ink)"}),
                                    dbc.Button(
                                        "Save Session",
                                        id="save-session-btn",
                                        color="success",
                                        className="mb-2 w-100",
                                    ),
                                    dcc.Upload(
                                        id="load-session-file",
                                        children=dbc.Button(
                                            "Load Session",
                                            color="secondary",
                                            outline=True,
                                            className="w-100",
                                        ),
                                        accept=".pkl",
                                    ),
                                ]
                            ),
                            className="card-accent-coral",
                        ),
                        width=4,
                    ),
                ],
                className="mb-3",
            ),
            dbc.Row(
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4(
                                    "👤 User Information (Required Before Running Analysis)",
                                    style={"color": "var(--ink)"},
                                ),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.Input(
                                                id="user-name",
                                                type="text",
                                                placeholder="Enter your name",
                                            ),
                                            width=6,
                                        ),
                                        dbc.Col(
                                            dbc.Input(
                                                id="user-id",
                                                type="text",
                                                placeholder="Enter user or lab ID",
                                            ),
                                            width=6,
                                        ),
                                    ]
                                ),
                                html.P(
                                    "This information is used to organize exported results "
                                    "into user-specific folders.",
                                    style={"fontSize": "13px", "color": "var(--muted)"},
                                    className="mt-2",
                                ),
                            ]
                        ),
                        className="card-accent-amber",
                    ),
                    width=12,
                ),
                className="mb-3",
            ),
            dbc.Row(dbc.Col(html.Div(id="home-status-message"), width=12), className="mb-3"),
            dbc.Row(
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4("📋 Quick Start Guide"),
                                html.Ol(
                                    [
                                        html.Li(
                                            "Select your analysis type (Proteomics, "
                                            "Metabolomics, or Phosphomatics)"
                                        ),
                                        html.Li(
                                            "Choose comparison method (T-Test or ANOVA multi-group)"
                                        ),
                                        html.Li("Upload your data file in the 'Data Upload' tab"),
                                        html.Li("Configure comparisons in the 'Comparisons' tab"),
                                        html.Li("Run analysis and visualize results"),
                                        html.Li("Perform enrichment analysis (optional)"),
                                        html.Li("Download all results in your preferred format"),
                                    ]
                                ),
                            ]
                        )
                    ),
                    width=12,
                )
            ),
        ]
    )


@callback(
    Output("store-app-config", "data"),
    Input("analysis-type", "value"),
    Input("comparison-method", "value"),
    Input("user-name", "value"),
    Input("user-id", "value"),
)
def sync_app_config(analysis_type, comparison_method, user_name, user_id):
    return {
        "analysis_type": analysis_type,
        "comparison_method": comparison_method,
        "user_name": user_name,
        "user_id": user_id,
    }


@callback(
    Output("home-status-message", "children"),
    Input("analysis-type", "value"),
    Input("comparison-method", "value"),
    Input("user-name", "value"),
    Input("user-id", "value"),
)
def render_status_message(analysis_type, comparison_method, user_name, user_id):
    user_name = (user_name or "").strip()
    user_id = (user_id or "").strip()
    has_user_info = bool(user_name) and bool(user_id)

    if analysis_type and comparison_method and has_user_info:
        return dbc.Alert(
            [
                html.Strong("Ready to proceed! "),
                f"Analysis: {analysis_type.upper()} | "
                f"Method: {_METHOD_LABELS.get(comparison_method, comparison_method)}",
                html.Br(),
                f"User: {user_name} | ID: {user_id}",
                html.Br(),
                dbc.Button(
                    "Go to Data Upload →",
                    id="goto-upload-btn",
                    color="primary",
                    size="sm",
                    className="mt-2",
                ),
            ],
            color="success",
        )
    if analysis_type and comparison_method and not has_user_info:
        return dbc.Alert(
            [
                html.Strong("User information required. "),
                "Please enter User Name and User ID before starting analysis.",
            ],
            color="warning",
        )
    return None
