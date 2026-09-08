import dash_bootstrap_components as dbc
import dash_daq as daq
import plotly.graph_objects as go
from dash import ALL, MATCH, Input, Output, State, callback, dcc, html
from dash.exceptions import PreventUpdate

from omics_app.data.columns import extract_group_names_from_columns

DEFAULT_PALETTE = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c"]


def layout() -> html.Div:
    return dbc.Row(
        [
            dbc.Col(
                dbc.Card(
                    dbc.CardBody(
                        [
                            html.H4("Assign Colors to Groups"),
                            html.P(
                                "Select colors for each experimental group. These colors "
                                "will be used consistently across all visualizations.",
                                style={"color": "var(--muted)", "fontSize": "14px"},
                            ),
                            html.Div(id="color-mapping-ui-container"),
                            dbc.Button(
                                "🎨 Save Color Mapping",
                                id="save-colors-btn",
                                color="success",
                                className="w-100 mt-3",
                            ),
                            html.Div(id="color-save-status", className="mt-2"),
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
                            html.H4("Color Preview"),
                            dcc.Graph(id="color-preview-plot", style={"height": "400px"}),
                            html.Hr(),
                            html.H5("Current Color Mapping:"),
                            html.Div(id="color-mapping-summary"),
                        ]
                    ),
                    className="card-accent-teal",
                ),
                width=6,
            ),
        ]
    )


def _sanitize(group: str) -> str:

    return "".join(c if c.isalnum() else "_" for c in group)


@callback(
    Output("color-mapping-ui-container", "children"),
    Input("store-main-data", "data"),
)
def populate_color_mapping_ui(main_data):

    if not main_data:
        raise PreventUpdate

    group_cols = main_data["column_index"]["group_cols"]
    groups = extract_group_names_from_columns(group_cols)

    if not groups:
        return html.P(
            style={"color": "var(--muted)"},
        )

    rows = []
    for i, group in enumerate(groups):
        default_color = DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)]
        swatch_id = {"type": "color-swatch-btn", "group": group}
        picker_wrapper_id = {"type": "color-picker-wrapper", "group": group}
        picker_id = {"type": "group-color-input", "group": group}
        hex_label_id = {"type": "color-hex-label", "group": group}

        rows.append(
            html.Div(
                [
                    html.Div(
                        [
                            html.Button(
                                id=swatch_id,
                                n_clicks=0,
                                style={
                                    "width": "32px",
                                    "height": "32px",
                                    "backgroundColor": default_color,
                                    "border": "2px solid var(--border)",
                                    "borderRadius": "8px",
                                    "cursor": "pointer",
                                    "padding": "0",
                                    "flexShrink": "0",
                                },
                            ),
                            html.Span(
                                group,
                                style={
                                    "fontWeight": "500",
                                    "marginLeft": "12px",
                                    "flex": "1",
                                },
                            ),
                            html.Span(
                                default_color,
                                id=hex_label_id,
                                style={
                                    "fontFamily": "monospace",
                                    "fontSize": "13px",
                                    "color": "var(--muted)",
                                },
                            ),
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": "4px"},
                    ),
                    # Picker is a plain div toggled via style (display
                    # none/block) rather than a dbc.Popover -- Popover's
                    # `target` matching against a pattern-matching dict ID
                    # is a code path I can't fully verify in this dbc
                    # version without a real browser; conditional display
                    # is the same mechanism already proven to work
                    # reliably for tab switching in app.py.
                    html.Div(
                        daq.ColorPicker(
                            id=picker_id,
                            value={"hex": default_color},
                            label="Pick a color",
                            style={"width": "220px"},
                        ),
                        id=picker_wrapper_id,
                        style={"display": "none", "marginTop": "8px"},
                    ),
                ],
                className="mb-2 p-2",
                style={
                    "border": "1px solid var(--border)",
                    "borderRadius": "8px",
                    "backgroundColor": "var(--surface)",
                },
            )
        )
    return rows


@callback(
    Output({"type": "color-picker-wrapper", "group": MATCH}, "style"),
    Input({"type": "color-swatch-btn", "group": MATCH}, "n_clicks"),
    State({"type": "color-picker-wrapper", "group": MATCH}, "style"),
    prevent_initial_call=True,
)
def toggle_color_picker(n_clicks, current_style):

    is_visible = (current_style or {}).get("display") == "block"
    new_style = dict(current_style or {})
    new_style["display"] = "none" if is_visible else "block"
    return new_style


@callback(
    Output({"type": "color-swatch-btn", "group": MATCH}, "style"),
    Output({"type": "color-hex-label", "group": MATCH}, "children"),
    Input({"type": "group-color-input", "group": MATCH}, "value"),
    State({"type": "color-swatch-btn", "group": MATCH}, "style"),
    prevent_initial_call=True,
)
def sync_swatch_with_picker(value, current_style):

    hex_color = (value or {}).get("hex", "#CCCCCC")
    new_style = dict(current_style or {})
    new_style["backgroundColor"] = hex_color
    return new_style, hex_color


@callback(
    Output("store-color-mapping", "data"),
    Output("color-save-status", "children"),
    Input("save-colors-btn", "n_clicks"),
    State({"type": "group-color-input", "group": ALL}, "value"),
    State({"type": "group-color-input", "group": ALL}, "id"),
    prevent_initial_call=True,
)
def save_colors(n_clicks, colors, color_ids):
    """Port of observeEvent(input$save_colors_btn, ...) (R line 2661)."""
    if not color_ids:
        return None, dbc.Alert("No groups available for color mapping yet.", color="warning")

    # daq.ColorPicker's value is {"hex": "#RRGGBB", "rgb": {...}} -- pull
    # just the hex string, matching R's plain hex-string color_mapping.
    color_mapping = {
        color_ids[i]["group"]: (colors[i] or {}).get("hex", DEFAULT_PALETTE[0])
        for i in range(len(color_ids))
    }
    return color_mapping, dbc.Alert("Color mapping saved!", color="success", className="py-2 mb-0")


@callback(
    Output("color-preview-plot", "figure"),
    Output("color-mapping-summary", "children"),
    Input("store-color-mapping", "data"),
)
def render_color_preview(color_mapping):
    """Port of output$color_preview_plot (R line 2681) and
    output$color_mapping_summary (R line 2704)."""
    if not color_mapping:
        empty_fig = go.Figure()
        empty_fig.update_layout(
            xaxis={"visible": False},
            yaxis={"visible": False},
            annotations=[
                {
                    "text": "Save a color mapping to see the preview",
                    "showarrow": False,
                    "font": {"color": "#9CA3AF"},
                }
            ],
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
        )
        return empty_fig, None

    groups = list(color_mapping.keys())
    colors = list(color_mapping.values())

    fig = go.Figure(go.Bar(x=groups, y=[1] * len(groups), marker_color=colors, showlegend=False))
    fig.update_layout(
        title="Group Color Assignment",
        xaxis={"tickangle": -45},
        yaxis={"visible": False},
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(t=40, b=80, l=20, r=20),
    )

    summary = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        style={
                            "width": "16px",
                            "height": "16px",
                            "backgroundColor": color,
                            "borderRadius": "3px",
                            "display": "inline-block",
                            "marginRight": "8px",
                            "border": "1px solid var(--border)",
                        }
                    ),
                    html.Span(f"{group}: ", style={"fontWeight": "600"}),
                    html.Span(color, style={"color": "var(--muted)", "fontFamily": "monospace"}),
                ],
                className="mb-2",
            )
            for group, color in color_mapping.items()
        ]
    )
    return fig, summary
