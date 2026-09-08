"""
Port of R tabPanel("📊 Visualization", ...) (app_12-02.R, lines
1132-1380) and its server logic (output$plot_volcano etc., R lines
5645-5886). See plotting/volcano.py, heatmap.py, pca.py, boxplot.py,
correlation.py -- the actual plot-building functions live there
(UI-free, testable) and this module is a thin Dash wrapper around
them, matching this project's established split (compare
stats/pipeline.py + ui/analysis.py).

Architectural note: R recomputes the full "results" object (intensity
matrix, significant_scaled, design, etc.) fresh for every plot from
rv$dea_results. This port instead reconstructs what each plot needs
directly from the flat records already stored in store-dea-results by
ui/analysis.py's run_analysis callback -- test_col_names/
control_col_names (added there specifically for this) give the
normalized-abundance columns that serve as R's intensity_matrix, and
Significance_pvalue/qvalue + sig_column give the "which rows are
significant" mask R gets from results$significant_pvalue.

Known gaps vs. R (see the individual plotting/*.py module docstrings
for the full list of simplifications within each plot):
  - Only the DEFAULT option of each R settings group is exposed here:
    volcano label modes are all wired, but heatmap "Top up + down"
    style variants, PCA's other 4 plot types (scree/loading/cumvar/
    variable), boxplot's raincloud/jitter styles, and 9 of
    correlation's 10 viz_types aren't in the dropdowns below yet --
    each is a real feature to add, not a bug, but they're not here.
  - No "Refresh All Plots" button (R line 1310) -- Dash callbacks
    already re-fire automatically on input change, so there's nothing
    for a manual refresh to do here.
  - ANOVA comparisons show a placeholder message in every plot slot,
    same as ui/analysis.py -- consistent with the rest of the port
    (see stats/anova_path.py's docstring for that broader gap).
"""

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, callback, dcc, html
from dash.exceptions import PreventUpdate

from omics_app.plotting.boxplot import build_boxplot_figure
from omics_app.plotting.correlation import build_correlation_figure
from omics_app.plotting.heatmap import build_heatmap_figure
from omics_app.plotting.pca import build_pca_figure
from omics_app.plotting.volcano import build_volcano_figure


def _empty_figure(message: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font=dict(size=14))
    fig.update_layout(xaxis_visible=False, yaxis_visible=False, template="plotly_white")
    return fig


def _load_comparison(dea_results: dict | None, comp_name: str | None):
    """
    Reconstructs (df, sample_columns, sample_groups, sig_mask,
    test_group_name, control_group_name) for one comparison from the
    flat dict stored by ui/analysis.py's run_analysis callback.
    Returns None if the comparison is missing, or isn't the 2-group
    ("normal") method this port supports.
    """
    if not dea_results or not comp_name or comp_name not in dea_results:
        return None
    entry = dea_results[comp_name]
    if "test_col_names" not in entry or "control_col_names" not in entry:
        # Either an ANOVA/uploaded-only entry, or from before this
        # tab's pipeline.py/analysis.py additions -- re-run Analysis
        # for this comparison to populate them.
        return None

    df = pd.DataFrame(entry["results"])
    rowname_col = entry.get("rowname_col")
    if rowname_col and rowname_col in df.columns:
        df = df.set_index(rowname_col)

    test_cols = entry["test_col_names"]
    control_cols = entry["control_col_names"]
    sample_groups = {c: entry["test_group_name"] for c in test_cols}
    sample_groups.update({c: entry["control_group_name"] for c in control_cols})

    sig_col = entry.get("sig_column", "Significance_pvalue")
    sig_mask = df[sig_col] if sig_col in df.columns else pd.Series(False, index=df.index)

    return {
        "df": df,
        "sample_columns": test_cols + control_cols,
        "sample_groups": sample_groups,
        "sig_mask": sig_mask,
        "test_group_name": entry["test_group_name"],
        "control_group_name": entry["control_group_name"],
    }


def layout() -> html.Div:
    settings = dbc.Card(
        dbc.CardBody(
            [
                html.H4("Plot Settings"),
                dcc.Dropdown(id="viz-comp-select", options=[], value=None, clearable=False),
                html.Hr(),
                html.H5("🌋 Volcano Plot"),
                dbc.Label("P-value cutoff:"),
                dbc.Input(
                    id="volcano-pval", type="number", value=0.05, min=0.001, max=0.1, step=0.005
                ),
                dbc.Label("Fold Change cutoff:", className="mt-2"),
                dbc.Input(id="volcano-fc", type="number", value=0, min=0, max=5, step=0.5),
                dbc.Label("Top genes to label:", className="mt-2"),
                dbc.Input(id="volcano-top-n", type="number", value=10, min=0, max=50),
                dbc.Label("Label Mode:", className="mt-2"),
                dcc.Dropdown(
                    id="volcano-label-mode",
                    options=[
                        {"label": "Protein symbol/name", "value": "symbol"},
                        {"label": "Accession IDs", "value": "accession"},
                        {"label": "Full protein names", "value": "fullname"},
                        {"label": "No labels", "value": "none"},
                    ],
                    value="fullname",
                    clearable=False,
                ),
                html.Hr(),
                html.H5("🔥 Heatmap"),
                dbc.Label("Heatmap Type:"),
                dcc.Dropdown(
                    id="heatmap-mode",
                    options=[
                        {"label": "All significant proteins", "value": "all"},
                        {"label": "Top upregulated proteins", "value": "up"},
                        {"label": "Top downregulated proteins", "value": "down"},
                        {"label": "Top up + down proteins", "value": "both"},
                    ],
                    value="all",
                    clearable=False,
                ),
                dbc.Label("Top proteins to display:", className="mt-2"),
                dbc.Input(id="heatmap-top-n", type="number", value=50, min=1, max=10000),
                dbc.Checklist(
                    id="heatmap-toggles",
                    options=[
                        {"label": "Scale data", "value": "scale"},
                        {"label": "Cluster rows (proteins)", "value": "cluster_rows"},
                        {"label": "Cluster columns (samples)", "value": "cluster_cols"},
                        {"label": "Show protein names", "value": "show_names"},
                    ],
                    value=["scale", "cluster_rows", "cluster_cols", "show_names"],
                    switch=True,
                    className="mt-2",
                ),
                dbc.Label("Row clustering distance:", className="mt-2"),
                dcc.Dropdown(
                    id="heatmap-row-distance",
                    options=[
                        {"label": d.capitalize(), "value": d}
                        for d in [
                            "euclidean",
                            "pearson",
                            "spearman",
                            "kendall",
                            "manhattan",
                            "canberra",
                            "maximum",
                        ]
                    ],
                    value="euclidean",
                    clearable=False,
                ),
                dbc.Label("Column clustering distance:", className="mt-2"),
                dcc.Dropdown(
                    id="heatmap-col-distance",
                    options=[
                        {"label": d.capitalize(), "value": d}
                        for d in [
                            "euclidean",
                            "pearson",
                            "spearman",
                            "kendall",
                            "manhattan",
                            "canberra",
                            "maximum",
                        ]
                    ],
                    value="euclidean",
                    clearable=False,
                ),
                html.Hr(),
                html.H5("📊 PCA Plot"),
                dbc.Label("PC X-axis:"),
                dbc.Input(id="pca-pc1", type="number", value=1, min=1, max=10),
                dbc.Label("PC Y-axis:", className="mt-2"),
                dbc.Input(id="pca-pc2", type="number", value=2, min=1, max=10),
                dbc.Checklist(
                    id="pca-toggles",
                    options=[
                        {"label": "Show 95% confidence ellipses", "value": "ellipses"},
                        {"label": "Show sample labels", "value": "labels"},
                    ],
                    value=["ellipses", "labels"],
                    switch=True,
                    className="mt-2",
                ),
                html.Hr(),
                html.H5("📦 Boxplot"),
                dbc.Label("Plot Style:"),
                dcc.Dropdown(
                    id="boxplot-style",
                    options=[
                        {"label": "Violin Plot", "value": "violin"},
                        {"label": "Traditional Boxplot", "value": "traditional"},
                    ],
                    value="violin",
                    clearable=False,
                ),
                html.Hr(),
                html.H5("🔗 Correlation"),
                dbc.Label("Method:"),
                dcc.Dropdown(
                    id="correlation-method",
                    options=[
                        {"label": "Pearson", "value": "pearson"},
                        {"label": "Spearman", "value": "spearman"},
                        {"label": "Kendall", "value": "kendall"},
                    ],
                    value="pearson",
                    clearable=False,
                ),
                dbc.Checklist(
                    id="correlation-toggles",
                    options=[{"label": "Show correlation coefficients", "value": "show_values"}],
                    value=["show_values"],
                    switch=True,
                    className="mt-2",
                ),
            ]
        ),
        className="card-accent-indigo",
    )

    def _plot_card(title, graph_id):
        placeholder = _empty_figure("Run Analysis on a comparison first.")
        return dbc.Card(
            dbc.CardBody([html.H5(title), dcc.Loading(dcc.Graph(id=graph_id, figure=placeholder))]),
            className="mb-3",
        )

    return dbc.Row(
        [
            dbc.Col(settings, width=3),
            dbc.Col(
                [
                    dbc.Row(
                        [
                            dbc.Col(_plot_card("Volcano Plot", "plot-volcano"), width=6),
                            dbc.Col(
                                _plot_card("Heatmap (Significant Proteins)", "plot-heatmap"),
                                width=6,
                            ),
                        ]
                    ),
                    dbc.Row(
                        [
                            dbc.Col(_plot_card("PCA Plot", "plot-pca"), width=6),
                            dbc.Col(
                                _plot_card("Boxplot (Intensity Distribution)", "plot-boxplot"),
                                width=6,
                            ),
                        ]
                    ),
                    dbc.Row(
                        dbc.Col(_plot_card("Correlation Heatmap", "plot-correlation"), width=12)
                    ),
                ],
                width=9,
            ),
        ]
    )


@callback(
    Output("viz-comp-select", "options"),
    Output("viz-comp-select", "value"),
    Input("store-dea-results", "data"),
)
def populate_viz_comparison_selector(dea_results):
    """Port of output$viz_comparison_selector (R line 5469)."""
    if not dea_results:
        return [], None
    options = [{"label": name, "value": name} for name in dea_results]
    return options, next(iter(dea_results))


@callback(
    Output("plot-volcano", "figure"),
    Input("viz-comp-select", "value"),
    Input("volcano-pval", "value"),
    Input("volcano-fc", "value"),
    Input("volcano-top-n", "value"),
    Input("volcano-label-mode", "value"),
    State("store-dea-results", "data"),
)
def render_volcano(comp_name, pval, fc, top_n, label_mode, dea_results):
    """Port of output$plot_volcano (R line 5645)."""
    loaded = _load_comparison(dea_results, comp_name)
    if loaded is None:
        raise PreventUpdate
    return build_volcano_figure(
        loaded["df"],
        pvalue_threshold=pval or 0.05,
        fc_threshold=fc or 0,
        top_n=top_n or 0,
        label_mode=label_mode or "fullname",
        test_group_name=loaded["test_group_name"],
        control_group_name=loaded["control_group_name"],
    )


@callback(
    Output("plot-heatmap", "figure"),
    Input("viz-comp-select", "value"),
    Input("heatmap-mode", "value"),
    Input("heatmap-top-n", "value"),
    Input("heatmap-toggles", "value"),
    Input("heatmap-row-distance", "value"),
    Input("heatmap-col-distance", "value"),
    State("store-dea-results", "data"),
    State("store-color-mapping", "data"),
)
def render_heatmap(
    comp_name, mode, top_n, toggles, row_distance, col_distance, dea_results, color_mapping
):
    """Port of output$plot_heatmap (R line 5680)."""
    loaded = _load_comparison(dea_results, comp_name)
    if loaded is None:
        raise PreventUpdate
    toggles = toggles or []
    return build_heatmap_figure(
        loaded["df"],
        loaded["sig_mask"],
        loaded["sample_columns"],
        loaded["sample_groups"],
        heatmap_mode=mode or "all",
        top_n=top_n or 50,
        scale_data="scale" in toggles,
        cluster_rows="cluster_rows" in toggles,
        cluster_cols="cluster_cols" in toggles,
        row_distance=row_distance or "euclidean",
        col_distance=col_distance or "euclidean",
        show_row_names="show_names" in toggles,
        color_mapping=color_mapping,
    )


@callback(
    Output("plot-pca", "figure"),
    Input("viz-comp-select", "value"),
    Input("pca-pc1", "value"),
    Input("pca-pc2", "value"),
    Input("pca-toggles", "value"),
    State("store-dea-results", "data"),
    State("store-color-mapping", "data"),
)
def render_pca(comp_name, pc1, pc2, toggles, dea_results, color_mapping):
    """Port of output$plot_pca (R line 5737), biplot branch only."""
    loaded = _load_comparison(dea_results, comp_name)
    if loaded is None:
        raise PreventUpdate
    toggles = toggles or []
    intensity_matrix = loaded["df"][loaded["sample_columns"]]
    return build_pca_figure(
        intensity_matrix,
        loaded["sample_groups"],
        pc1=int(pc1 or 1),
        pc2=int(pc2 or 2),
        show_ellipses="ellipses" in toggles,
        show_labels="labels" in toggles,
        color_mapping=color_mapping,
    )


@callback(
    Output("plot-boxplot", "figure"),
    Input("viz-comp-select", "value"),
    Input("boxplot-style", "value"),
    State("store-dea-results", "data"),
    State("store-color-mapping", "data"),
)
def render_boxplot(comp_name, style, dea_results, color_mapping):
    """Port of output$plot_boxplot (R line 5812)."""
    loaded = _load_comparison(dea_results, comp_name)
    if loaded is None:
        raise PreventUpdate
    intensity_matrix = loaded["df"][loaded["sample_columns"]]
    return build_boxplot_figure(
        intensity_matrix,
        loaded["sample_groups"],
        plot_style=style or "violin",
        color_mapping=color_mapping,
    )


@callback(
    Output("plot-correlation", "figure"),
    Input("viz-comp-select", "value"),
    Input("correlation-method", "value"),
    Input("correlation-toggles", "value"),
    State("store-dea-results", "data"),
)
def render_correlation(comp_name, method, toggles, dea_results):
    """Port of output$plot_correlation (R line 5851)."""
    loaded = _load_comparison(dea_results, comp_name)
    if loaded is None:
        raise PreventUpdate
    intensity_matrix = loaded["df"][loaded["sample_columns"]]
    return build_correlation_figure(
        intensity_matrix,
        method=method or "pearson",
        show_values="show_values" in (toggles or []),
    )
