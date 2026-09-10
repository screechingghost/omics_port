"""
Dash app shell -- mirrors the 9-tab navbarPage structure of app_12-02.R
(ui <- navbarPage(...), line 450), restructured as a sidebar-nav layout.
"""

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html

from omics_app.server import cache
from omics_app.ui import (
    analysis,
    colors,
    comparisons,
    download,
    enrichment,
    home,
    session,
    upload,
    visualization,
)

app = dash.Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.FLATLY,
        "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap",
    ],
    suppress_callback_exceptions=True,
)
app.title = "Multi-Omics Analysis"
cache.init_app(app.server)

# --- Design tokens ---------------------------------------------------------
# bg:      #F7F8FC  main content background
# surface: #FFFFFF  card background
# ink:     #1A1F36  primary text
# muted:   #6B7280  secondary text
# indigo:  #4F46E5  brand / primary actions / active nav item
# teal:    #0EA5A0  secondary accent
# coral:   #F0654B  tertiary accent
# amber:   #E8A23D  quaternary accent
# border:  #E4E7F0  hairline dividers
#
# Structural change from the previous pass: color now lives in ONE
# place (the sidebar) rather than being split between a gradient hero
# AND a gradient tab bar -- two loud elements were competing. The
# per-card accent borders on the Home tab stay, since those encode
# real grouping information rather than pure decoration.
app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            :root {
                --bg: #F7F8FC;
                --surface: #FFFFFF;
                --ink: #1A1F36;
                --muted: #6B7280;
                --indigo: #4F46E5;
                --indigo-dark: #241E63;
                --teal: #0EA5A0;
                --coral: #F0654B;
                --amber: #E8A23D;
                --border: #E4E7F0;
            }
            html, body { height: 100%; margin: 0; }
            body {
                background-color: var(--bg) !important;
                color: var(--ink);
                font-family: 'Inter', -apple-system, sans-serif;
            }
            h1, h2, h3, h4, h5, h6 {
                font-family: 'Inter', sans-serif;
                font-weight: 600;
                letter-spacing: -0.01em;
                color: var(--ink);
            }
            p, label, .form-label { color: var(--ink); }
            small, .text-muted, .form-text { color: var(--muted) !important; }

            /* --- App shell: fixed sidebar + scrollable content --- */
            .app-shell { display: flex; height: 100vh; overflow: hidden; }
            .sidebar {
                width: 248px;
                height: 100vh;
                flex-shrink: 0;
                overflow-y: auto;
                background: linear-gradient(180deg, var(--indigo) 0%, var(--indigo-dark) 100%);
                padding: 20px 0;
                display: flex;
                flex-direction: column;
            }
            .sidebar-brand {
                padding: 4px 20px 20px 24px;
                border-bottom: 1px solid rgba(255,255,255,0.15);
                margin-bottom: 12px;
            }
            .sidebar-brand-title {
                color: #FFFFFF;
                font-size: 16px;
                font-weight: 600;
                letter-spacing: -0.01em;
                margin: 0;
            }
            .sidebar-brand-subtitle {
                color: rgba(255,255,255,0.65);
                font-size: 12px;
                margin: 2px 0 0 0;
            }
            .main-content {
                flex: 1;
                padding: 32px 40px;
                height: 100vh;
                overflow-y: auto;
                min-width: 0;
                max-width: 1200px;
            }

            /* --- Vertical nav (dbc.Tabs vertical=True) --- */
            #main-tabs .nav-link {
                border: none;
                border-left: 3px solid transparent;
                border-radius: 0;
                color: rgba(255,255,255,0.72);
                font-weight: 500;
                font-size: 14px;
                text-align: left;
                padding: 11px 20px 11px 21px;
                background: transparent;
                transition: background 0.15s ease, color 0.15s ease;
            }
            #main-tabs .nav-link:hover {
                color: #FFFFFF;
                background: rgba(255,255,255,0.06);
            }
            #main-tabs .nav-link.active {
                color: #FFFFFF;
                font-weight: 600;
                background: rgba(255,255,255,0.10);
                border-left-color: var(--amber);
            }

            /* --- Cards --- */
            .card {
                border: 1px solid var(--border);
                border-radius: 10px;
                box-shadow: 0 1px 2px rgba(26, 31, 54, 0.04);
            }
            .card-accent-indigo { border-top: 3px solid var(--indigo); }
            .card-accent-teal   { border-top: 3px solid var(--teal); }
            .card-accent-coral  { border-top: 3px solid var(--coral); }
            .card-accent-amber  { border-top: 3px solid var(--amber); }

            .btn-primary {
                background-color: var(--indigo);
                border-color: var(--indigo);
                font-weight: 500;
            }
            .btn-primary:hover {
                background-color: #4038C7;
                border-color: #4038C7;
            }
            .btn-success {
                background-color: var(--teal);
                border-color: var(--teal);
                font-weight: 500;
            }
            .btn-success:hover {
                background-color: #0C8983;
                border-color: #0C8983;
            }

            .form-control {
                border-color: var(--border) !important;
                color: var(--ink) !important;
            }
            .form-control:focus {
                border-color: var(--indigo) !important;
                box-shadow: 0 0 0 3px rgba(79, 70, 229, 0.12) !important;
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""

TABS = [
    ("home", "🏠 Home"),
    ("upload", "📁 Data Upload"),
    ("comparisons", "⚖️ Comparisons"),
    ("colors", "🎨 Color Mapping"),
    ("analysis", "🔬 Analysis"),
    ("visualization", "📊 Visualization"),
    ("enrichment", "🧬 Enrichment"),
    ("download", "💾 Download"),
    ("session", "ℹ️ Session Info"),
]

# tab_id -> layout function. Add entries here as each tab gets built.
TAB_LAYOUTS = {
    "home": home.layout,
    "upload": upload.layout,
    "comparisons": comparisons.layout,
    "colors": colors.layout,
    "analysis": analysis.layout,
    "visualization": visualization.layout,
    "enrichment": enrichment.layout,
    "download": download.layout,
    "session": session.layout,
}

app.layout = html.Div(
    [
        # rv-equivalent: server-side application state.
        # If main_data/dea_results get large, back this with flask-caching
        # (server-side session) rather than the default client-side JSON
        # store -- see README "State management" section.
        dcc.Store(id="store-main-data"),
        dcc.Store(id="store-app-config"),
        dcc.Store(id="store-comparisons"),
        dcc.Store(id="store-dea-results"),
        dcc.Store(id="store-enrichment-results"),
        dcc.Store(id="store-color-mapping"),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.P("Multi-Omics Suite", className="sidebar-brand-title"),
                                html.P(
                                    "Proteomics · Metabolomics · Phosphomatics",
                                    className="sidebar-brand-subtitle",
                                ),
                            ],
                            className="sidebar-brand",
                        ),
                        dbc.Tabs(
                            [dbc.Tab(label=label, tab_id=tab_id) for tab_id, label in TABS],
                            id="main-tabs",
                            active_tab="home",
                            className="flex-column",
                        ),
                    ],
                    className="sidebar",
                ),
                # Every tab's layout is built once, up front, and stays in
                # the DOM permanently -- switching tabs only toggles which
                # panel is visible (see set_active_panel below). This is
                # what fixes state loss on tab switch: dropdown selections,
                # loaded previews, etc. used to be destroyed and rebuilt
                # from scratch every time render_tab() replaced tab-content
                # wholesale, even though the underlying dcc.Store data was
                # never actually lost.
                html.Div(
                    [
                        html.Div(
                            layout_fn(),
                            id=f"tab-panel-{tab_id}",
                            style={"display": "block" if tab_id == "home" else "none"},
                        )
                        for tab_id, layout_fn in TAB_LAYOUTS.items()
                    ]
                    + [
                        html.Div(
                            f"TODO: build layout for tab '{tab_id}'",
                            id=f"tab-panel-{tab_id}",
                            style={"display": "none"},
                        )
                        for tab_id, _label in TABS
                        if tab_id not in TAB_LAYOUTS
                    ],
                    className="main-content",
                ),
            ],
            className="app-shell",
        ),
    ]
)


@dash.callback(
    [Output(f"tab-panel-{tab_id}", "style") for tab_id, _label in TABS],
    Input("main-tabs", "active_tab"),
)
def set_active_panel(active_tab: str):
    """Shows the active tab's panel, hides the rest. Panels themselves
    are never destroyed (see app.layout above), so any state a user has
    already set -- dropdown picks, a loaded data preview, filled-in
    comparisons -- survives switching away and back."""
    return [
        {"display": "block"} if tab_id == active_tab else {"display": "none"}
        for tab_id, _label in TABS
    ]


# "Go to Data Upload →" button on the Home tab jumps tabs, mirroring
# observeEvent(input$goto_upload, ...) (R line 1892).
@dash.callback(
    Output("main-tabs", "active_tab"),
    Input("goto-upload-btn", "n_clicks"),
    prevent_initial_call=True,
)
def goto_upload(n_clicks):
    return "upload"


if __name__ == "__main__":
    app.run(debug=True)
