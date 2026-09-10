"""
Port of R tabPanel("ℹ️ Session Info", ...) (app_12-02.R, lines 1630-1659).

UI-ONLY, as requested -- these three panels are static placeholders.
R populates them from session_info_text/analysis_log_text/
packages_info_text renderers (base R sessionInfo(), an in-memory log
of user actions, and loaded-package versions); none of that exists on
the Python side yet, so each panel just says so instead of showing
fabricated data.
"""

import dash_bootstrap_components as dbc
from dash import html


def _info_card(title: str, placeholder: str) -> dbc.Card:
    return dbc.Card(
        dbc.CardBody(
            [
                html.H3(title, style={"color": "var(--ink)"}),
                html.Pre(
                    placeholder,
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
    return html.Div(
        [
            _info_card("Session Information", "Not available yet -- this will report Python/package versions once wired up."),
            _info_card("Analysis Log", "No actions logged yet -- this will show a running log of uploads, comparisons, and analyses run this session."),
            _info_card("Loaded Packages", "Not available yet -- this will list the installed package versions in use."),
        ]
    )
