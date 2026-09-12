# pyright: reportCallIssue=false, reportInvalidTypeForm=false
import dash_bootstrap_components as dbc
import pandas as pd
from dash import Input, Output, State, callback, dash_table, dcc, html
from dash.exceptions import PreventUpdate

from omics_app.stats.pipeline import AnalysisError, run_two_group_analysis


def layout() -> html.Div:
    return dbc.Row(
        [
            dbc.Col(
                dbc.Card(
                    dbc.CardBody(
                        [
                            html.H4("Select Comparison for Analysis"),
                            dcc.Dropdown(id="analysis-comp-select", options=[], value=None),
                            html.Hr(),
                            html.H5("📊 Filtering Parameters"),
                            html.Div(
                                [
                                    dbc.Label("Minimum Valid Values (%):"),
                                    dcc.Slider(
                                        id="min-valid-percent",
                                        min=30,
                                        max=100,
                                        step=5,
                                        value=70,
                                        marks={i: str(i) for i in range(30, 101, 10)},
                                    ),
                                    html.Small(
                                        "Lower = more proteins kept but lower quality. "
                                        "Recommended: 50-70% discovery, 70-80% validation.",
                                        className="text-muted",
                                    ),
                                ],
                                className="p-2 mb-2",
                                style={"backgroundColor": "#FFF3CD", "borderRadius": "6px"},
                            ),
                            html.H5("📈 Significance Thresholds"),
                            html.Div(
                                [
                                    dbc.Label("P-value cutoff:"),
                                    dbc.Input(
                                        id="pvalue-threshold",
                                        type="number",
                                        value=0.05,
                                        min=0.001,
                                        max=0.1,
                                        step=0.01,
                                    ),
                                    dbc.Label("Use for filtering:", className="mt-2"),
                                    dcc.Dropdown(
                                        id="significance-method",
                                        options=[
                                            {"label": "FDR-adjusted (q-value)", "value": "fdr"},
                                            {"label": "Raw p-value", "value": "raw"},
                                        ],
                                        value="raw",
                                        clearable=False,
                                    ),
                                    html.Hr(),
                                    html.H6("🔄 Fold Change Filter (Optional)"),
                                    dbc.Label("Log2 Fold Change threshold:"),
                                    dbc.Input(
                                        id="fc-threshold",
                                        type="number",
                                        value=0,
                                        min=0,
                                        max=5,
                                        step=0.1,
                                    ),
                                    dbc.Label("Apply filter:", className="mt-2"),
                                    dcc.Dropdown(
                                        id="fc-operator",
                                        options=[
                                            {"label": "No fold change filter", "value": "none"},
                                            {"label": "Greater than or equal (≥)", "value": "gte"},
                                            {"label": "Greater than (>)", "value": "gt"},
                                            {"label": "Less than or equal (≤)", "value": "lte"},
                                            {"label": "Less than (<)", "value": "lt"},
                                            {"label": "Absolute value ≥", "value": "abs"},
                                        ],
                                        value="none",
                                        clearable=False,
                                    ),
                                ],
                                className="p-2 mb-2",
                                style={"backgroundColor": "#E3F2FD", "borderRadius": "6px"},
                            ),
                            html.Hr(),
                            dbc.Button(
                                "▶ Run Analysis",
                                id="run-analysis-btn",
                                color="primary",
                                size="lg",
                                className="w-100 mb-2",
                            ),
                            dbc.Button(
                                "⏩ Run All Comparisons",
                                id="run-all-analysis-btn",
                                color="success",
                                className="w-100",
                            ),
                            html.Hr(),
                            html.H5("Analysis Progress:"),
                            dcc.Loading(html.Div(id="analysis-progress")),
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
                            html.H4("Analysis Results Summary"),
                            html.Div(id="analysis-results-summary"),
                            html.Hr(),
                            html.H5("Significant Results Preview:"),
                            dcc.Loading(html.Div(id="significant-results-table")),
                        ]
                    ),
                    className="card-accent-teal",
                ),
                width=8,
            ),
        ]
    )


@callback(
    Output("analysis-comp-select", "options"),
    Output("analysis-comp-select", "value"),
    Input("store-comparisons", "data"),
)
def populate_analysis_comparison_selector(comparisons_data):

    comparisons = (comparisons_data or {}).get("comparisons") or []
    if not comparisons:
        return [], None
    options = [{"label": c["name"], "value": c["name"]} for c in comparisons]
    return options, comparisons[0]["name"]


def _build_top5_table(result: dict):

    df = result["results_df"]
    pvalue_sig_df = df[df["Significance_pvalue"]]
    if pvalue_sig_df.empty:
        return None

    available_cols = [c for c in ["logFC", "P.Value", "adj.P.Val"] if c in pvalue_sig_df.columns]
    if not available_cols:
        return html.P("Column information not available", style={"color": "var(--muted)"})

    top5 = pvalue_sig_df.sort_values("P.Value").head(5)[available_cols].reset_index()
    return dash_table.DataTable(
        data=top5.to_dict("records"),
        columns=[{"name": c, "id": c} for c in top5.columns],
        style_table={"overflowX": "auto"},
        style_cell={
            "fontSize": "12px",
            "fontFamily": "monospace",
            "textAlign": "left",
            "padding": "6px",
        },
        style_header={"textAlign": "center", "fontWeight": "600"},
    )


def _build_summary(result: dict, comp_name: str, params: dict) -> html.Div:

    fc_op_labels = {"gte": "≥", "gt": ">", "lte": "≤", "lt": "<", "abs": "|FC| ≥"}
    fc_line = "Fold Change Filter: None (all fold changes included)"
    if params["fc_operator"] != "none" and params["fc_threshold"] > 0:
        op_label = fc_op_labels.get(params["fc_operator"], params["fc_operator"])
        fc_line = f"Fold Change Filter: log2FC {op_label} {params['fc_threshold']:.2f}"

    sig_method_label = (
        "FDR-adjusted (q-value)" if params["significance_method"] == "fdr" else "Raw p-value"
    )

    stat_blocks = [
        ("Proteins Analyzed", result["n_kept"], "indigo"),
        ("Significant (p-value)", result["n_significant_pvalue"], "teal"),
        ("Significant (q-value)", result["n_significant_qvalue"], "coral"),
        ("Upregulated", result["n_upregulated"], "amber"),
        ("Downregulated", result["n_downregulated"], "muted"),
    ]

    top5_table = _build_top5_table(result)
    top5_section = (
        [html.Hr(), html.P("Top 5 Most Significant", className="mb-1 fw-bold"), top5_table]
        if top5_table is not None
        else []
    )

    return html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        html.Div(
                            [
                                html.Div(
                                    str(value),
                                    style={
                                        "fontSize": "24px",
                                        "fontWeight": "700",
                                        "color": f"var(--{accent})"
                                        if accent != "muted"
                                        else "var(--muted)",
                                    },
                                ),
                                html.Div(
                                    label, style={"fontSize": "11px", "color": "var(--muted)"}
                                ),
                            ],
                            style={"textAlign": "center"},
                        ),
                        width=True,
                    )
                    for label, value, accent in stat_blocks
                ],
                className="mb-2",
            ),
            html.Hr(),
            html.P([html.Strong("Comparison: "), comp_name], className="mb-1"),
            html.P(
                [
                    html.Strong("Groups: "),
                    result["test_group_name"],
                    " vs ",
                    result["control_group_name"],
                ],
                className="mb-1",
            ),
            html.Hr(),
            html.P("Analysis Parameters", className="mb-1 fw-bold"),
            html.P(f"Min Valid Values: {params['min_valid_percent']}%", className="mb-1"),
            html.P(f"P-value Threshold: {params['pvalue_threshold']:.3f}", className="mb-1"),
            html.P(f"Significance Method: {sig_method_label}", className="mb-1"),
            html.P(fc_line, className="mb-1"),
            *top5_section,
        ]
    )


def _build_results_table(result: dict):
    df = result["results_df"]
    sig_df = df[df[result["sig_column"]]]
    if sig_df.empty:
        return dbc.Alert("No proteins reached the significance threshold.", color="warning")

    display_cols = [
        c for c in ["logFC", "P.Value", "adj.P.Val", "Comparison"] if c in sig_df.columns
    ]
    preview = sig_df[display_cols].reset_index().head(100)
    return dash_table.DataTable(
        data=preview.to_dict("records"),
        columns=[{"name": c, "id": c} for c in preview.columns],
        page_size=10,
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_cell={
            "fontSize": "12px",
            "fontFamily": "monospace",
            "textAlign": "left",
            "padding": "8px",
        },
        style_header={"textAlign": "center", "fontWeight": "600"},
        sort_action="native",
    )


@callback(
    Output("store-dea-results", "data", allow_duplicate=True),
    Output("analysis-progress", "children"),
    Output("analysis-results-summary", "children"),
    Output("significant-results-table", "children"),
    Input("run-analysis-btn", "n_clicks"),
    State("analysis-comp-select", "value"),
    State("store-comparisons", "data"),
    State("store-main-data", "data"),
    State("store-app-config", "data"),
    State("min-valid-percent", "value"),
    State("pvalue-threshold", "value"),
    State("significance-method", "value"),
    State("fc-threshold", "value"),
    State("fc-operator", "value"),
    State("store-dea-results", "data"),
    prevent_initial_call=True,
)
def run_analysis(
    n_clicks,
    comp_name,
    comparisons_data,
    main_data,
    app_config,
    min_valid_percent,
    pvalue_threshold,
    significance_method,
    fc_threshold,
    fc_operator,
    existing_dea_results,
):

    # User-info gate, matching R lines 4342-4349
    user_name = ((app_config or {}).get("user_name") or "").strip()
    user_id = ((app_config or {}).get("user_id") or "").strip()
    if not user_name or not user_id:
        alert = dbc.Alert(
            "Please enter User Name and User ID on the Home tab before running analysis.",
            color="danger",
        )
        return existing_dea_results, alert, None, None

    if not comp_name or not main_data:
        raise PreventUpdate

    comparisons = (comparisons_data or {}).get("comparisons") or []
    comparison = next((c for c in comparisons if c["name"] == comp_name), None)
    if comparison is None:
        raise PreventUpdate

    if comparison.get("method") != "normal":
        alert = dbc.Alert(
            "ANOVA analysis isn't wired up yet in this port -- only 2-group "
            "('normal') comparisons can be run here so far.",
            color="warning",
        )
        return existing_dea_results, alert, None, None

    df = pd.DataFrame(main_data["data"])
    rowname_col = main_data.get("rowname_col")
    if rowname_col and rowname_col in df.columns:
        df = df.set_index(rowname_col)

    try:
        result = run_two_group_analysis(
            df,
            comparison,
            min_valid_percent=min_valid_percent or 70,
            pvalue_threshold=pvalue_threshold or 0.05,
            significance_method=significance_method or "raw",
            fc_threshold=fc_threshold or 0,
            fc_operator=fc_operator or "none",
        )
    except AnalysisError as e:
        alert = dbc.Alert(str(e), color="danger")
        return existing_dea_results, alert, None, None

    dea_results = dict(existing_dea_results or {})
    dea_results[comp_name] = {
        "results": result["results_df"].reset_index().to_dict("records"),
        "n_kept": result["n_kept"],
        "n_significant": result["n_significant"],
        "test_group_name": result["test_group_name"],
        "control_group_name": result["control_group_name"],
        "sig_column": result["sig_column"],
        # For ui/visualization.py: which normalized-abundance columns
        # belong to which group, and the index column name needed to
        # restore `results` back into a properly-indexed DataFrame
        # (reset_index() above names that column after rowname_col).
        "test_col_names": result["test_col_names"],
        "control_col_names": result["control_col_names"],
        "rowname_col": rowname_col,
    }

    progress = dbc.Alert(
        f"Analysis complete for '{comp_name}'.", color="success", className="py-2 mb-0"
    )
    params = {
        "min_valid_percent": min_valid_percent or 70,
        "pvalue_threshold": pvalue_threshold or 0.05,
        "significance_method": significance_method or "raw",
        "fc_threshold": fc_threshold or 0,
        "fc_operator": fc_operator or "none",
    }
    summary = _build_summary(result, comp_name, params)
    table = _build_results_table(result)
    return dea_results, progress, summary, table
