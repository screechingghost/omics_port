"""
Port of R tabPanel("⚖️ Comparisons", ...) (app_12-02.R, lines 968-1010)
and its server logic (metadata_cols_selector: 2053, comparisons_setup_ui:
2222, setup_comparisons_btn handler: 2324, auto-detect abundance columns
from group selection: 2480-2542).

Known simplification vs. the R app: the "no 'Found in Sample Group'
columns detected, enter everything manually" fallback path (R lines
2271-2287, 2306-2316) isn't ported yet -- this assumes group_cols is
non-empty, which is the common case for real Proteome Discoverer
exports. Add the manual-entry fallback if you hit real data without
those columns.
"""

import re

import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, callback, dash_table, dcc, html
from dash.exceptions import PreventUpdate

from omics_app.data.columns import extract_group_names_from_columns


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
    """
    IMPORTANT: every dropdown below explicitly sets value=None (or []
    for multi-selects), even though that's the default anyway. This
    card gets rebuilt from scratch whenever store-main-data changes
    (a new file is uploaded) -- if a dropdown's `value` prop were
    simply omitted instead of explicitly reset, Dash treats "prop not
    present in this render" as "don't change it", not "reset it". A
    dropdown that previously had a value selected from an older
    upload would then keep that stale value while receiving a NEW
    `options` list (from the new file) that may no longer contain it.
    react-select then fails to find the matching option object to
    display and throws "Cannot read properties of null (reading
    'label')" in the browser. Explicit value=None/[] here closes that
    gap regardless of whether Dash/React happens to reuse the
    underlying component instance across renders.
    """
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
    else:  # "normal" 2-group path (default)
        body = html.Div(
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
# Port of the observeEvent pair at R lines 2505-2535: when a group is
# picked, pre-fill the abundance-column selector by matching the group
# name as a substring of the abundance column names.


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
    Input("setup-comparisons-btn", "n_clicks"),
    State({"type": "comp-name", "index": ALL}, "value"),
    State({"type": "comp-name", "index": ALL}, "id"),
    State({"type": "comp-group1", "index": ALL}, "value"),
    State({"type": "comp-group2", "index": ALL}, "value"),
    State({"type": "comp-abundance-g1", "index": ALL}, "value"),
    State({"type": "comp-abundance-g2", "index": ALL}, "value"),
    State({"type": "comp-anova-groups", "index": ALL}, "value"),
    State("metadata-cols-checklist", "value"),
    State("store-app-config", "data"),
    State("store-main-data", "data"),
    prevent_initial_call=True,
)
def setup_comparisons(
    n_clicks,
    names,
    name_ids,
    group1s,
    group2s,
    abund_g1s,
    abund_g2s,
    anova_groups,
    metadata_cols,
    app_config,
    main_data,
):
    """Port of observeEvent(input$setup_comparisons_btn, ...) (R line 2324).

    Also precomputes each comparison's selected_cols (R line 2380-2397:
    metadata + group + abundance columns actually used), matching the R
    app's approach of freezing rv$comparison_list at Setup time rather
    than recomputing it live when previewed.

    Known simplification for the ANOVA path: R auto-derives per-group
    abundance columns from "Found in Sample Group" flags (not yet
    ported -- see comparisons_setup_ui docstring); the preview for
    ANOVA comparisons here shows metadata + the selected group columns
    only, not their underlying abundance columns.
    """
    method = (app_config or {}).get("comparison_method") or "normal"
    metadata_cols = metadata_cols or []
    all_columns = list(main_data["data"][0].keys()) if main_data and main_data.get("data") else []

    comparisons = []
    for idx, name in enumerate(names):
        comp_index = name_ids[idx]["index"]
        entry = {"index": comp_index, "name": name, "method": method}

        if method == "anova":
            groups = anova_groups[idx] if idx < len(anova_groups) else []
            entry["groups"] = groups
            wanted_cols = metadata_cols + (groups or [])
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
    return payload, "\n".join(summary_lines)


# --- Comparison Preview: pick a comparison, see the column-sliced data
# that would actually feed into it. Port of output$comparison_preview_selector
# (R line 2624) and output$comparison_preview_table (R line 2630).


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
