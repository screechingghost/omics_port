"""
Port of the enrichment plot renderer (app_12-02.R, generate_enrichment_plot
and friends, ~line 10250 onward) -- all 9 of R's plot_type options are
built here, adapted to this port's Enrichr-based results shape
(Direction/Source/Term/Overlap/PValue/AdjPValue/Genes -- see
stats/enrichment.py) rather than R's clusterProfiler S4 result objects:

  - barplot, dotplot: direct ports (R lines 11077-11323).
  - go_ontology: R lines 10916-11074, grouping by ONTOLOGY (BP/MF/CC) --
    here grouped by Source instead (our 3 GO Enrichr libraries stand in
    for R's 3 ontology aspects), Up/Down side by side instead of
    ggarrange'd.
  - bubble: R lines 11324-11373, faceted by Category x Direction -- here
    faceted by Source x Direction using Plotly subplots.
  - up_down_compare: R lines 11375-11405, signed horizontal bar chart --
    direct port.
  - cnetplot, emapplot (R lines 11407-11454): R uses `enrichplot`'s
    igraph-based layouts on the actual clusterProfiler S4 objects. There's
    no equivalent Python package wired into this project (no networkx
    dependency either), so these are rebuilt from the flat results_df as
    simple manual-layout networks (circular gene-term bipartite graph for
    cnetplot, circular term-similarity graph via Jaccard overlap for
    emapplot) -- visually different from `enrichplot`'s force-directed
    layouts, but convey the same "which genes/terms are connected"
    relationship. Flagged as an approximation in each docstring.
  - chord_up, chord_down (R lines 11456+, via `circlize`): true circular
    chord diagrams have no Plotly equivalent. Rebuilt as a Sankey diagram
    (genes -> pathways) instead -- same many-to-many relationship, common
    substitute for a chord diagram, but visually a flow diagram rather
    than a circular one. Also flagged as an approximation.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

UP_COLOR = "#CC3311"
DOWN_COLOR = "#4477AA"
ONTOLOGY_COLORS = {
    "GO_Biological_Process": "#298C8C",
    "GO_Molecular_Function": "#A00000",
    "GO_Cellular_Component": "#B8B8B8",
    "Reactome_Pathways": "#7A5195",
}
ONTOLOGY_LABELS = {
    "GO_Biological_Process": "Biological Process",
    "GO_Molecular_Function": "Molecular Function",
    "GO_Cellular_Component": "Cellular Component",
    "Reactome_Pathways": "Reactome Pathway",
}


def _message_figure(message: str) -> go.Figure:
    """Port of R's enrichment_message_plot() -- shown when there's not
    enough data for a given plot type rather than erroring."""
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font={"size": 14})
    fig.update_layout(xaxis_visible=False, yaxis_visible=False, template="plotly_white")
    return fig


def _parse_overlap(overlap: str) -> tuple[int, int]:
    try:
        k, n = str(overlap).split("/")
        return int(k), int(n)
    except (ValueError, AttributeError):
        return 0, 1


def _truncate(text: str, max_len: int) -> str:
    if pd.isna(text):
        return ""

    text = str(text)
    return text if len(text) <= max_len else text[: max_len - 3] + "..."


def _split_genes(genes_field) -> list[str]:
    """Enrichr's Genes column is semicolon-separated (gseapy's own
    formatting, not R's comma-separated clusterProfiler geneID field)."""
    if not genes_field or pd.isna(genes_field):
        return []
    return [g.strip() for g in str(genes_field).split(";") if g.strip()]


def build_enrichment_barplot(df: pd.DataFrame, top_n: int = 15) -> go.Figure:
    if df is None or df.empty:
        return _message_figure("No enrichment results to display")

    data = df.sort_values("PValue").head(top_n).copy()
    data["neglog10p"] = -np.log10(data["PValue"].clip(lower=1e-300))
    data = data.iloc[::-1]  # smallest p-value ends up drawn at the top

    fig = go.Figure()
    for direction, color in [("Up", UP_COLOR), ("Down", DOWN_COLOR)]:
        subset = data[data["Direction"] == direction]
        if subset.empty:
            continue
        fig.add_trace(
            go.Bar(
                x=subset["neglog10p"],
                y=subset["Term"],
                orientation="h",
                name=direction,
                marker_color=color,
                hovertemplate="%{y}<br>-log10(p)=%{x:.2f}<extra></extra>",
            )
        )

    fig.update_layout(
        title=f"Top {min(top_n, len(df))} Enriched Terms",
        xaxis_title="-Log10(p-value)",
        barmode="group",
        template="plotly_white",
        margin={"l": 300, "t": 50},
    )
    return fig


def build_enrichment_dotplot(df: pd.DataFrame, top_n: int = 15) -> go.Figure:
    if df is None or df.empty:
        return _message_figure("No enrichment results to display")

    data = df.sort_values("PValue").head(top_n).copy()
    parsed = data["Overlap"].apply(_parse_overlap)
    data["k"] = [p[0] for p in parsed]
    data["n"] = [p[1] for p in parsed]
    data["GeneRatio"] = data["k"] / data["n"].replace(0, 1)
    data = data.iloc[::-1]

    max_k = max(int(data["k"].max()), 1)
    fig = go.Figure(
        go.Scatter(
            x=data["GeneRatio"],
            y=data["Term"],
            mode="markers",
            marker={
                "size": data["k"],
                "sizemode": "area",
                "sizeref": 2.0 * max_k / (40.0**2),
                "sizemin": 4,
                "color": data["PValue"],
                "colorscale": "Viridis_r",
                "showscale": True,
                "colorbar": {"title": "p-value"},
            },
            text=data["Direction"],
            hovertemplate="%{y}<br>GeneRatio=%{x:.2f}<br>%{text}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Top {min(top_n, len(df))} Enriched Terms",
        xaxis_title="Gene Ratio",
        template="plotly_white",
        margin={"l": 300, "t": 50},
    )
    return fig


def build_go_ontology_plot(df: pd.DataFrame, top_n_per_ontology: int = 10) -> go.Figure:
    """Port of generate_go_ontology_plot (R lines 10916-11074). Groups by
    Source (our stand-in for R's GO ONTOLOGY: BP/MF/CC) instead of a
    dedicated ONTOLOGY column, and renders Up/Down as side-by-side
    subplots instead of ggarrange -- same information, no gtable
    dependency."""
    go_df = df[df["Source"].isin(ONTOLOGY_LABELS) & (df["Source"] != "Reactome_Pathways")]
    if go_df.empty:
        return _message_figure("No significant GO enrichment results")

    directions = [d for d in ("Up", "Down") if not go_df[go_df["Direction"] == d].empty]
    if not directions:
        return _message_figure("No significant GO enrichment results")

    fig = make_subplots(
        rows=1,
        cols=len(directions),
        subplot_titles=[f"{d}-Regulated GO Terms" for d in directions],
    )
    for col, direction in enumerate(directions, start=1):
        subset = go_df[go_df["Direction"] == direction]
        parts = []
        for source in ONTOLOGY_LABELS:
            top = subset[subset["Source"] == source].sort_values("PValue").head(top_n_per_ontology)
            parts.append(top)
        combined = pd.concat(parts) if parts else subset.iloc[0:0]
        if combined.empty:
            continue
        combined = combined.assign(
            k=[_parse_overlap(o)[0] for o in combined["Overlap"]],
            OntologyLabel=combined["Source"].map(ONTOLOGY_LABELS),
            TermLabel=combined["Term"].map(lambda t: _truncate(t, 45)),
        ).sort_values(["OntologyLabel", "k"], ascending=[True, False])

        for source, label in ONTOLOGY_LABELS.items():
            group = combined[combined["Source"] == source]
            if group.empty:
                continue
            fig.add_trace(
                go.Bar(
                    x=group["TermLabel"],
                    y=group["k"],
                    name=label,
                    legendgroup=label,
                    showlegend=(col == 1),
                    marker_color=ONTOLOGY_COLORS[source],
                    text=group["k"],
                    textposition="outside",
                ),
                row=1,
                col=col,
            )

    fig.update_xaxes(tickangle=90)
    fig.update_layout(
        title="GO Terms by Ontology",
        yaxis_title="Gene Count",
        template="plotly_white",
        barmode="group",
        legend_title_text="Gene Ontology",
        margin={"b": 160},
    )
    return fig


def build_enrichment_bubble_plot(df: pd.DataFrame, max_per_type: int = 15) -> go.Figure:
    """Port of generate_enrichment_bubble_plot (R lines 11324-11373).
    R facets by Category x Direction with ggplot's facet_grid; built here
    as a Plotly subplot grid (rows = Source, cols = Direction)."""
    if df is None or df.empty:
        return _message_figure("No significant enrichment results for bubble plot")

    parts = [
        group.sort_values("PValue").head(max_per_type)
        for _, group in df.groupby(["Source", "Direction"])
    ]
    data = pd.concat(parts) if parts else df.iloc[0:0]
    if data.empty:
        return _message_figure("No significant enrichment results for bubble plot")

    parsed = data["Overlap"].apply(_parse_overlap)
    data = data.assign(
        k=[p[0] for p in parsed],
        n=[p[1] for p in parsed],
        TermLabel=data["Term"].map(lambda t: _truncate(t, 55)),
    )
    data["GeneRatio"] = data["k"] / data["n"].replace(0, 1)

    sources = [s for s in ONTOLOGY_LABELS if s in data["Source"].unique()]
    directions = [d for d in ("Up", "Down") if d in data["Direction"].unique()]
    if not sources or not directions:
        return _message_figure("No significant enrichment results for bubble plot")

    fig = make_subplots(
        rows=len(sources),
        cols=len(directions),
        column_titles=directions,
        row_titles=[ONTOLOGY_LABELS[s] for s in sources],
        shared_xaxes=False,
    )
    max_k = max(int(data["k"].max()), 1)
    for r, source in enumerate(sources, start=1):
        for c, direction in enumerate(directions, start=1):
            cell = data[(data["Source"] == source) & (data["Direction"] == direction)]
            if cell.empty:
                continue
            cell = cell.sort_values("GeneRatio")
            fig.add_trace(
                go.Scatter(
                    x=cell["GeneRatio"],
                    y=cell["TermLabel"],
                    mode="markers",
                    marker={
                        "size": cell["k"],
                        "sizemode": "area",
                        "sizeref": 2.0 * max_k / (30.0**2),
                        "sizemin": 4,
                        "color": cell["AdjPValue"],
                        "colorscale": "Viridis_r",
                        "showscale": (r == 1 and c == 1),
                        "colorbar": {"title": "Adj. p-value"},
                    },
                    showlegend=False,
                    hovertemplate="%{y}<br>GeneRatio=%{x:.2f}<br>Count=%{marker.size}<extra></extra>",
                ),
                row=r,
                col=c,
            )

    fig.update_layout(
        title="Enrichment Bubble Plot",
        template="plotly_white",
        height=max(420, 140 * len(sources)),
        margin={"l": 220},
    )
    return fig


def build_up_down_comparison_plot(df: pd.DataFrame, max_per_type: int = 12) -> go.Figure:
    """Port of generate_up_down_comparison_plot (R lines 11375-11405) --
    signed horizontal bar chart, Up bars to the right, Down bars to the
    left, sorted by signed gene count."""
    if df is None or df.empty:
        return _message_figure("No significant enrichment results for up vs down comparison")

    parts = [
        group.sort_values("PValue").head(max_per_type)
        for _, group in df.groupby(["Source", "Direction"])
    ]
    data = pd.concat(parts) if parts else df.iloc[0:0]
    if data.empty:
        return _message_figure("No significant enrichment results for up vs down comparison")

    parsed = data["Overlap"].apply(_parse_overlap)
    data = data.assign(k=[p[0] for p in parsed])
    data["SignedCount"] = np.where(data["Direction"] == "Down", -data["k"], data["k"])
    data["TermLabel"] = data.apply(
        lambda r: f"{_truncate(r['Term'], 55)} | {ONTOLOGY_LABELS.get(r['Source'], r['Source'])}",
        axis=1,
    )
    data = data.sort_values("SignedCount")

    fig = go.Figure()
    for direction, color in [("Up", UP_COLOR), ("Down", DOWN_COLOR)]:
        subset = data[data["Direction"] == direction]
        if subset.empty:
            continue
        fig.add_trace(
            go.Bar(
                x=subset["SignedCount"],
                y=subset["TermLabel"],
                orientation="h",
                name=direction,
                marker_color=color,
                hovertemplate="%{y}<br>Count=%{customdata}<extra></extra>",
                customdata=subset["k"],
            )
        )
    fig.add_vline(x=0, line_color="black", line_width=1)
    max_abs = max(int(data["k"].max()), 1)
    fig.update_layout(
        title="Up vs Down Regulated Enrichment Comparison",
        xaxis_title="Gene Count",
        xaxis={"tickvals": [-max_abs, 0, max_abs], "ticktext": [str(max_abs), "0", str(max_abs)]},
        legend_title_text="Regulation",
        template="plotly_white",
        margin={"l": 340},
        height=max(420, 26 * len(data)),
    )
    return fig


def _circular_positions(n: int, radius: float = 1.0, offset: float = 0.0):
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False) + offset
    return radius * np.cos(angles), radius * np.sin(angles)


def build_enrichment_cnetplot(df: pd.DataFrame, max_terms: int = 8) -> go.Figure:
    """Approximation of generate_enrichment_cnetplot (R lines 11407-11427),
    which calls enrichplot::cnetplot on a clusterProfiler S4 object with
    igraph's force-directed layout. No such package is wired into this
    port, so this rebuilds the same gene-term bipartite graph -- term
    nodes on an inner ring, their genes on an outer ring, edges linking
    each term to its overlap genes -- from the flat Genes column instead.
    Same "which genes drive which terms" relationship, different
    (deterministic, not physics-based) layout.

    R prefers Reactome terms when available (prefer_reactome=TRUE) --
    matched here too.
    """
    if df is None or df.empty:
        return _message_figure("No significant enrichment terms available for cnetplot")

    reactome = df[df["Source"] == "Reactome_Pathways"]
    pool = reactome if not reactome.empty else df
    terms = pool.sort_values("PValue").head(max_terms).reset_index(drop=True)
    if terms.empty:
        return _message_figure("No significant enrichment terms available for cnetplot")

    term_genes = {row["Term"]: _split_genes(row.get("Genes")) for _, row in terms.iterrows()}
    all_genes = sorted({g for genes in term_genes.values() for g in genes})
    if not all_genes:
        return _message_figure("No gene-level detail available for cnetplot")

    term_labels = [_truncate(t, 40) for t in term_genes]
    term_x, term_y = _circular_positions(len(term_labels), radius=1.0)
    gene_x, gene_y = _circular_positions(len(all_genes), radius=2.4)
    gene_pos = dict(zip(all_genes, zip(gene_x, gene_y)))

    edge_x, edge_y = [], []
    for i, (term, genes) in enumerate(term_genes.items()):
        for g in genes:
            gx, gy = gene_pos[g]
            edge_x += [term_x[i], gx, None]
            edge_y += [term_y[i], gy, None]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line={"color": "rgba(150,150,150,0.5)", "width": 1},
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=gene_x,
            y=gene_y,
            mode="markers+text",
            text=all_genes,
            textposition="top center",
            textfont={"size": 9},
            marker={"size": 8, "color": "#B8B8B8"},
            name="Gene",
            hovertemplate="%{text}<extra></extra>",
        )
    )
    term_sizes = [8 + 4 * len(g) for g in term_genes.values()]
    fig.add_trace(
        go.Scatter(
            x=term_x,
            y=term_y,
            mode="markers+text",
            text=term_labels,
            textposition="middle center",
            textfont={"size": 10, "color": "white"},
            marker={"size": term_sizes, "color": UP_COLOR},
            name="Term",
            hovertemplate="%{text}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Gene-Pathway Network ({'Reactome' if not reactome.empty else 'Top terms'})",
        template="plotly_white",
        showlegend=True,
        xaxis={"visible": False},
        yaxis={"visible": False},
        height=560,
    )
    return fig


def build_enrichment_map_plot(
    df: pd.DataFrame, max_terms: int = 20, similarity_threshold: float = 0.2
) -> go.Figure:
    """Approximation of generate_enrichment_map_plot (R lines 11429-11454),
    which runs enrichplot::pairwise_termsim + emapplot (igraph
    force-directed layout weighted by gene-set Jaccard similarity). Same
    idea here -- edges between terms sharing >= similarity_threshold of
    their genes (Jaccard index) -- but positioned on a circle rather than
    a force-directed layout, since no graph-layout package is wired into
    this project."""
    if df is None or df.empty:
        return _message_figure("No significant enrichment terms available for enrichment map")

    terms = df.sort_values("PValue").head(max_terms).reset_index(drop=True)
    if len(terms) < 2:
        return _message_figure("At least two significant terms are required for enrichment map")

    gene_sets = [set(_split_genes(g)) for g in terms["Genes"]]
    labels = [_truncate(t, 40) for t in terms["Term"]]
    parsed = terms["Overlap"].apply(_parse_overlap)
    counts = [p[0] for p in parsed]
    x, y = _circular_positions(len(terms))

    edge_x, edge_y = [], []
    for i in range(len(terms)):
        for j in range(i + 1, len(terms)):
            a, b = gene_sets[i], gene_sets[j]
            if not a or not b:
                continue
            jaccard = len(a & b) / len(a | b)
            if jaccard >= similarity_threshold:
                edge_x += [x[i], x[j], None]
                edge_y += [y[i], y[j], None]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line={"color": "rgba(150,150,150,0.6)", "width": 1},
            hoverinfo="skip",
            showlegend=False,
        )
    )
    max_count = max(counts) if counts else 1
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers+text",
            text=labels,
            textposition="top center",
            textfont={"size": 9},
            marker={
                "size": [10 + 20 * (c / max_count) for c in counts],
                "color": terms["AdjPValue"],
                "colorscale": "Viridis_r",
                "showscale": True,
                "colorbar": {"title": "Adj. p-value"},
            },
            hovertemplate="%{text}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Enrichment Map (term similarity)",
        template="plotly_white",
        xaxis={"visible": False},
        yaxis={"visible": False},
        height=560,
    )
    return fig


def build_enrichment_chord_plot(df: pd.DataFrame, direction: str, max_terms: int = 10) -> go.Figure:
    """Approximation of R's Reactome-specific circlize chord diagram
    (generate_enrichment_chord, R lines 11456+) for one regulation
    direction. A true chord/circos diagram (circular layout, ribbons
    between outer arc segments) has no Plotly primitive; this uses a
    Sankey diagram instead -- same many-to-many gene<->pathway
    relationship, drawn as left-to-right flow rather than a circle.
    Restricted to Reactome results, matching R's scope for this plot."""
    label = "Up-Regulated" if direction == "up" else "Down-Regulated"
    reactome = df[
        (df["Source"] == "Reactome_Pathways")
        & (df["Direction"] == ("Up" if direction == "up" else "Down"))
    ]
    if reactome.empty:
        return _message_figure(f"No significant Reactome pathways for {label.lower()} genes")

    terms = reactome.sort_values("PValue").head(max_terms).reset_index(drop=True)
    pathway_genes = {row["Term"]: _split_genes(row.get("Genes")) for _, row in terms.iterrows()}
    genes = sorted({g for gs in pathway_genes.values() for g in gs})
    if not genes:
        return _message_figure(f"No gene-level detail available for {label.lower()} pathways")

    pathway_labels = [_truncate(t, 40) for t in pathway_genes]
    node_labels = genes + pathway_labels
    gene_index = {g: i for i, g in enumerate(genes)}
    pathway_index = {t: len(genes) + i for i, t in enumerate(pathway_genes)}

    sources, targets, values = [], [], []
    for term, gs in pathway_genes.items():
        for g in gs:
            sources.append(gene_index[g])
            targets.append(pathway_index[term])
            values.append(1)

    color = UP_COLOR if direction == "up" else DOWN_COLOR
    fig = go.Figure(
        go.Sankey(
            node={
                "label": node_labels,
                "pad": 10,
                "thickness": 14,
                "color": ["#B8B8B8"] * len(genes) + [color] * len(pathway_labels),
            },
            link={
                "source": sources,
                "target": targets,
                "value": values,
                "color": "rgba(150,150,150,0.4)",
            },
        )
    )
    fig.update_layout(
        title=f"Pathway Chord Diagram (Sankey approximation) - {label}",
        template="plotly_white",
        height=max(420, 22 * len(node_labels)),
    )
    return fig
