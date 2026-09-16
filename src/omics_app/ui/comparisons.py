# pyright: reportCallIssue=false, reportInvalidTypeForm=false
import re

import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, callback, dash_table, dcc, html
from dash.exceptions import PreventUpdate

from omics_app.data.columns import extract_group_names_from_columns
from omics_app.stats.session_log import append_log_entry


def layout() -> html.Div:
    return dbc.Row(
        [
            dbc.Col(
                dbc.Card(
                    dbc.CardBody(
                        [
                            html.H4("Comparison Configuration"),
                            dbc.Label("Number of Comparisons:"),
                            dbc.Input(
                                id="n-comparisons",
                                type="number",
                                value=1,
                                min=1,
                                max=10,
                                step=1,
                            ),
                            html.Hr(),
                            html.H5("Metadata Columns (Common for all comparisons):"),
                            html.Div(
                                [
                                    dbc.Button(
                                        "Select All",
                                        id="metadata-select-all-btn",
                                        size="sm",
                                        color="success",
                                        outline=True,
                                        className="me-2",
                                    ),
                                    dbc.Button(
                                        "Clear",
                                        id="metadata-clear-btn",
                                        size="sm",
                                        color="secondary",
                                        outline=True,
                                    ),
                                ],
                                className="mb-2",
                            ),
                            dcc.Checklist(
                                id="metadata-cols-checklist",
                                options=[],
                                value=[],
                                labelStyle={"display": "block"},
                            ),
                            dbc.Button(
                                "Setup Comparisons",
                                id="setup-comparisons-btn",
                                color="primary",
                                className="w-100 mt-3",
                            ),
                        ]
                    )
                ),
                width=4,
            ),
            dbc.Col(
                [
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4("Comparison Details"),
                                html.Div(id="comparisons-setup-container"),
                            ]
                        ),
                        className="mb-3",
                    ),
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4("Comparison Summary"),
                                html.Pre(id="comparison-summary-text"),
                            ]
                        ),
                        className="mb-3",
                    ),
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H4("Comparison Preview"),
                                dbc.Label("Select Comparison to Preview:"),
                                dcc.Dropdown(
                                    id="comparison-preview-select",
                                    options=[],
                                    value=None,
                                    clearable=False,
                                ),
                                dcc.Loading(
                                    html.Div(id="comparison-preview-table", className="mt-3")
                                ),
                            ]
                        ),
                        className="card-accent-teal",
                    ),
                ],
                width=8,
            ),
        ]
    )


# --- Metadata column checklist -------------------------------------------


@callback(
    Output("metadata-cols-checklist", "options"),
    Output("metadata-cols-checklist", "value"),
    Input("store-main-data", "data"),
)
def populate_metadata_checklist(main_data):
    if not main_data:
        raise PreventUpdate
    candidates = main_data["column_index"]["metadata_candidates"]
    return [{"label": c, "value": c} for c in candidates], candidates


@callback(
    Output("metadata-cols-checklist", "value", allow_duplicate=True),
    Input("metadata-select-all-btn", "n_clicks"),
    State("store-main-data", "data"),
    prevent_initial_call=True,
)
def select_all_metadata(n_clicks, main_data):
    if not main_data:
        raise PreventUpdate
    return main_data["column_index"]["metadata_candidates"]


@callback(
    Output("metadata-cols-checklist", "value", allow_duplicate=True),
    Input("metadata-clear-btn", "n_clicks"),
    prevent_initial_call=True,
)
def clear_metadata(n_clicks):
    return []


# --- Dynamic per-comparison cards -----------------------------------------


def _comparison_card(i: int, group_cols: list[str], abundance_cols: list[str], method: str):

    header = html.Div(
        [
            html.H5(f"Comparison {i}", style={"color": "#2c3e50"}),
            dbc.Label("Comparison Name:"),
            dbc.Input(
                id={"type": "comp-name", "index": i},
                type="text",
                value=f"Comparison_{i}",
            ),
        ]
    )

    if method == "anova":
        if group_cols:
            body = html.Div(
                [
                    dbc.Label("Select Groups:"),
                    dcc.Dropdown(
                        id={"type": "comp-anova-groups", "index": i},
                        options=[{"label": g, "value": g} for g in group_cols],
                        value=[],  # explicit reset -- see note above _comparison_card
                        multi=True,
                    ),
                ]
            )
        else:
            # Port of R lines 2253-2260: no "Found in Sample Group" columns
            # detected. R offers a full manual per-group builder here too
            # (numericInput + dynamic group rows); that manual builder
            # still isn't wired in this port -- stats/pipeline.py's
            # run_multi_group_analysis (see README) always auto-matches
            # abundance columns from selected group names, so without
            # group columns to select from, there's nothing to feed it.
            body = html.Div(
                dbc.Alert(
                    [
                        html.I(className="fas fa-info-circle me-1"),
                        (
                            "No 'Found in Sample Group' columns detected. Manual per-group "
                            "abundance assignment for ANOVA isn't available yet in this port -- "
                            "see README Known Gaps."
                        ),
                    ],
                    color="info",
                    className="py-2 mb-0",
                    style={"fontSize": "13px"},
                )
            )
    else:  # "normal" 2-group path (default)
        if group_cols:
            group_inputs = html.Div(
                [
                    dbc.Label("Group 1 (Test):"),
                    dcc.Dropdown(
                        id={"type": "comp-group1", "index": i},
                        options=[{"label": g, "value": g} for g in group_cols],
                        value=None,  # explicit reset -- see note above _comparison_card
                    ),
                    dbc.Label("Group 2 (Control):", className="mt-2"),
                    dcc.Dropdown(
                        id={"type": "comp-group2", "index": i},
                        options=[{"label": g, "value": g} for g in group_cols],
                        value=None,
                    ),
                ]
            )
        else:
            # Manual fallback -- port of R lines 2271-2287. No "Found in
            # Sample Group" columns exist, so the group NAME can't be
            # auto-detected from a column header; the user types plain
            # names instead. Deliberately separate component "type" from
            # comp-group1/comp-group2 (not just a differently-populated
            # version of the same id) so the auto_fill_abundance_g1/g2
            # callbacks below -- which assume a group value is a column
            # header to pattern-match against -- never fire for manual
            # entries; R's own manual path has no such auto-detection
            # either, abundance columns are always picked explicitly here.
            #
            # No pipeline.py/stats changes needed for this: a manually
            # typed name flows into comparison["group1"/"group2"] exactly
            # like an auto-detected column header does --
            # extract_group_names_from_columns() already returns a plain
            # string unchanged when it doesn't match the "Found in
            # Sample...Group...:" pattern (data/columns.py).
            group_inputs = html.Div(
                [
                    dbc.Alert(
                        [
                            html.I(className="fas fa-info-circle me-1"),
                            "No 'Found in Sample Group' columns detected. Create Test and Control "
                            "group names manually here for color mapping, analysis, and plots.",
                        ],
                        color="info",
                        className="py-2 mb-2",
                        style={"fontSize": "13px"},
                    ),
                    dbc.Label("Group 1 Name (Test):"),
                    dbc.Input(
                        id={"type": "comp-group1-manual", "index": i},
                        type="text",
                        value="Test",
                    ),
                    dbc.Label("Group 2 Name (Control):", className="mt-2"),
                    dbc.Input(
                        id={"type": "comp-group2-manual", "index": i},
                        type="text",
                        value="Control",
                    ),
                ]
            )

        body = html.Div(
            [
                group_inputs,
                dbc.Label("Abundance Columns (Group 1 - Test):", className="mt-2"),
                dcc.Dropdown(
                    id={"type": "comp-abundance-g1", "index": i},
                    options=[{"label": c, "value": c} for c in abundance_cols],
                    value=[],
                    multi=True,
                ),
                dbc.Label("Abundance Columns (Group 2 - Control):", className="mt-2"),
                dcc.Dropdown(
                    id={"type": "comp-abundance-g2", "index": i},
                    options=[{"label": c, "value": c} for c in abundance_cols],
                    value=[],
                    multi=True,
                ),
            ]
        )

    return dbc.Card(
        dbc.CardBody([header, html.Hr(), body]),
        className="mb-3",
        style={"backgroundColor": "#f0f8ff"},
    )


@callback(
    Output("comparisons-setup-container", "children"),
    Input("n-comparisons", "value"),
    Input("store-main-data", "data"),
    Input("store-app-config", "data"),
)
def render_comparison_cards(n_comparisons, main_data, app_config):
    """Port of output$comparisons_setup_ui (R line 2222)."""
    if not main_data or not n_comparisons or n_comparisons < 1:
        raise PreventUpdate

    column_index = main_data["column_index"]
    group_cols = column_index["group_cols"]
    abundance_cols = column_index["abundance_cols"]
    method = (app_config or {}).get("comparison_method") or "normal"

    if not abundance_cols:
        return dbc.Alert(
            [
                html.Strong("No abundance columns found! "),
                "Looking for columns matching: 'Abundance:'",
                html.Br(),
                "Your columns: ",
                html.Small(", ".join(list(main_data["data"][0].keys())[:20])),
            ],
            color="warning",
        )

    return [
        _comparison_card(i, group_cols, abundance_cols, method)
        for i in range(1, int(n_comparisons) + 1)
    ]


# --- Auto-detect abundance columns from group selection --------------------


def _auto_match_abundance(group_col: str | None, abundance_cols: list[str]) -> list[str]:
    if not group_col:
        raise PreventUpdate
    group_name = extract_group_names_from_columns([group_col])[0]
    if not group_name:
        raise PreventUpdate
    pattern = re.escape(group_name)
    return [c for c in abundance_cols if re.search(pattern, c, re.IGNORECASE)]


@callback(
    Output({"type": "comp-abundance-g1", "index": ALL}, "value"),
    Input({"type": "comp-group1", "index": ALL}, "value"),
    State("store-main-data", "data"),
    prevent_initial_call=True,
)
def auto_fill_abundance_g1(group1_values, main_data):
    if not main_data:
        raise PreventUpdate
    abundance_cols = main_data["column_index"]["abundance_cols"]
    return [_auto_match_abundance(g, abundance_cols) if g else [] for g in group1_values]


@callback(
    Output({"type": "comp-abundance-g2", "index": ALL}, "value"),
    Input({"type": "comp-group2", "index": ALL}, "value"),
    State("store-main-data", "data"),
    prevent_initial_call=True,
)
def auto_fill_abundance_g2(group2_values, main_data):
    if not main_data:
        raise PreventUpdate
    abundance_cols = main_data["column_index"]["abundance_cols"]
    return [_auto_match_abundance(g, abundance_cols) if g else [] for g in group2_values]


# --- Setup Comparisons button: collect everything into store-comparisons ---


@callback(
    Output("store-comparisons", "data", allow_duplicate=True),
    Output("comparison-summary-text", "children"),
    Output("store-analysis-log", "data", allow_duplicate=True),
    Input("setup-comparisons-btn", "n_clicks"),
    State({"type": "comp-name", "index": ALL}, "value"),
    State({"type": "comp-name", "index": ALL}, "id"),
    State({"type": "comp-group1", "index": ALL}, "value"),
    State({"type": "comp-group2", "index": ALL}, "value"),
    State({"type": "comp-group1-manual", "index": ALL}, "value"),
    State({"type": "comp-group2-manual", "index": ALL}, "value"),
    State({"type": "comp-abundance-g1", "index": ALL}, "value"),
    State({"type": "comp-abundance-g2", "index": ALL}, "value"),
    State({"type": "comp-anova-groups", "index": ALL}, "value"),
    State("metadata-cols-checklist", "value"),
    State("store-app-config", "data"),
    State("store-main-data", "data"),
    State("store-analysis-log", "data"),
    prevent_initial_call=True,
)
def setup_comparisons(
    n_clicks,
    names,
    name_ids,
    group1s,
    group2s,
    group1_manuals,
    group2_manuals,
    abund_g1s,
    abund_g2s,
    anova_groups,
    metadata_cols,
    app_config,
    main_data,
    existing_log,
):

    method = (app_config or {}).get("comparison_method") or "normal"
    metadata_cols = metadata_cols or []
    all_columns = list(main_data["data"][0].keys()) if main_data and main_data.get("data") else []
    group_cols_available = bool(((main_data or {}).get("column_index") or {}).get("group_cols"))

    comparisons = []
    skipped: list[str] = []
    for idx, name in enumerate(names):
        comp_index = name_ids[idx]["index"]
        entry = {"index": comp_index, "name": name, "method": method}

        if method == "anova":
            groups = anova_groups[idx] if idx < len(anova_groups) else []
            entry["groups"] = groups
            wanted_cols = metadata_cols + (groups or [])
        elif not group_cols_available:
            # Manual 2-group fallback -- port of R lines 2352-2387,
            # including its exact validation rules (required names,
            # required abundance selections on both sides, no column
            # double-assigned to both groups).
            group1_manual = (
                (group1_manuals[idx] if idx < len(group1_manuals) else "") or ""
            ).strip()
            group2_manual = (
                (group2_manuals[idx] if idx < len(group2_manuals) else "") or ""
            ).strip()
            abundance_g1 = abund_g1s[idx] if idx < len(abund_g1s) else []
            abundance_g2 = abund_g2s[idx] if idx < len(abund_g2s) else []
            abundance_g1, abundance_g2 = abundance_g1 or [], abundance_g2 or []

            if not group1_manual or not group2_manual:
                skipped.append(f"{name}: skipped -- Test and Control group names are required.")
                continue
            if not abundance_g1 or not abundance_g2:
                skipped.append(
                    f"{name}: skipped -- select abundance columns for both Test and Control."
                )
                continue
            if set(abundance_g1) & set(abundance_g2):
                skipped.append(
                    f"{name}: skipped -- some abundance columns are assigned to both Test and Control."
                )
                continue

            entry["group1"] = group1_manual
            entry["group2"] = group2_manual
            entry["abundance_g1"] = abundance_g1
            entry["abundance_g2"] = abundance_g2
            wanted_cols = metadata_cols + abundance_g1 + abundance_g2
        else:
            group1 = group1s[idx] if idx < len(group1s) else None
            group2 = group2s[idx] if idx < len(group2s) else None
            abundance_g1 = abund_g1s[idx] if idx < len(abund_g1s) else []
            abundance_g2 = abund_g2s[idx] if idx < len(abund_g2s) else []
            entry["group1"] = group1
            entry["group2"] = group2
            entry["abundance_g1"] = abundance_g1 or []
            entry["abundance_g2"] = abundance_g2 or []
            wanted_cols = (
                metadata_cols + [group1, group2] + (abundance_g1 or []) + (abundance_g2 or [])
            )

        # unique() preserving order, filtered to columns that actually
        # exist in main_data -- mirrors R's selected_cols[selected_cols
        # %in% colnames(rv$main_data)]
        seen = set()
        selected_cols = []
        for col in wanted_cols:
            if col and col in all_columns and col not in seen:
                seen.add(col)
                selected_cols.append(col)
        entry["selected_cols"] = selected_cols

        comparisons.append(entry)

    payload = {"method": method, "metadata_cols": metadata_cols, "comparisons": comparisons}

    summary_lines = [f"{len(comparisons)} comparison(s) configured ({method}):"]
    for c in comparisons:
        if method == "anova":
            summary_lines.append(f"  {c['name']}: groups = {c.get('groups')}")
        else:
            summary_lines.append(
                f"  {c['name']}: {c.get('group1')} (n={len(c.get('abundance_g1') or [])}) "
                f"vs {c.get('group2')} (n={len(c.get('abundance_g2') or [])})"
            )
    if skipped:
        summary_lines.append("")
        summary_lines.append("Skipped:")
        summary_lines.extend(f"  {line}" for line in skipped)

    return (
        payload,
        "\n".join(summary_lines),
        append_log_entry(
            existing_log, f"Comparisons configured: {len(comparisons)} comparison(s) ({method})"
        ),
    )


# --- Comparison Preview: pick a comparison, see the column-sliced data


@callback(
    Output("comparison-preview-select", "options"),
    Output("comparison-preview-select", "value"),
    Input("store-comparisons", "data"),
)
def populate_preview_selector(comparisons_data):
    """Port of output$comparison_preview_selector (R line 2624)."""
    comparisons = (comparisons_data or {}).get("comparisons") or []
    if not comparisons:
        return [], None
    options = [{"label": c["name"], "value": c["name"]} for c in comparisons]
    return options, comparisons[0]["name"]


@callback(
    Output("comparison-preview-table", "children"),
    Input("comparison-preview-select", "value"),
    State("store-comparisons", "data"),
    State("store-main-data", "data"),
)
def render_comparison_preview(selected_name, comparisons_data, main_data):
    """
    Port of output$comparison_preview_table (R line 2630): shows the
    first 50 rows of main_data restricted to that comparison's
    selected_cols (metadata + group + abundance columns actually used) --
    lets the user sanity-check the comparison before running DEA on it.
    """
    if not selected_name or not comparisons_data or not main_data:
        raise PreventUpdate

    comparisons = comparisons_data.get("comparisons") or []
    entry = next((c for c in comparisons if c["name"] == selected_name), None)
    if entry is None:
        raise PreventUpdate

    selected_cols = entry.get("selected_cols") or []
    if not selected_cols:
        return dbc.Alert(
            "No columns resolved for this comparison -- check group/abundance "
            "selections and click Setup Comparisons again.",
            color="warning",
        )

    rows = main_data["data"][:50]
    preview_rows = [{col: row.get(col) for col in selected_cols} for row in rows]

    return dash_table.DataTable(
        data=preview_rows,
        columns=[{"name": c, "id": c} for c in selected_cols],
        page_size=10,
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_cell={"fontSize": "12px", "fontFamily": "monospace", "textAlign": "left"},
        style_header={"textAlign": "center", "fontWeight": "600"},
    )
