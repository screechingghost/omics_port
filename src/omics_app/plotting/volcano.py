"""
Port of generate_volcano_plot_advanced (app_12-02.R, line 6418) --
2-group path only. The ANOVA branch (R lines 6427-6486, a plain
F-statistic vs. -log10(p) plot for multi-group data) isn't ported --
consistent with the rest of this codebase, which doesn't have a
working ANOVA pipeline yet (see stats/anova_path.py, ui/analysis.py).

Known simplifications vs. R:
- R renders with EnhancedVolcano (Bioconductor), which draws connector
  lines from labels to points and auto-avoids overlaps. This builds a
  Plotly scatter with the same color/cutoff logic and static
  annotation labels for the top N up/down proteins by logFC --
  visually close, not pixel-identical, and labels can overlap at high
  top_n on crowded plots (no repel algorithm).
- R's coloring test (line 6560-6562, strict '<') and its Expression
  classification used for top-gene selection (line 6526-6528, '<=')
  use slightly different pvalue comparisons, so a handful of R's
  points sitting exactly on pval_cutoff can be labeled but colored
  gray. This port uses one consistent '<=' test everywhere -- a
  deliberate cleanup, not a carried-over quirk.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

UP_COLOR = "#BB0A1E"
DOWN_COLOR = "royalblue"
UNCHANGED_COLOR = "#AAAAAA"


def _make_unique(values: pd.Series) -> pd.Series:
    """Mirrors R's make.unique(..., sep = "_")."""
    seen: dict[str, int] = {}
    out = []
    for v in values:
        if v not in seen:
            seen[v] = 0
            out.append(v)
        else:
            seen[v] += 1
            out.append(f"{v}_{seen[v]}")
    return pd.Series(out, index=values.index)


def _pick_labels(df: pd.DataFrame, label_mode: str) -> pd.Series:
    """Port of the label_values selection (R lines 6490-6515)."""
    if "Accession" in df.columns:
        accession = df["Accession"].astype(str)
    else:
        accession = pd.Series(df.index.astype(str), index=df.index)

    if "Description" in df.columns:
        description = df["Description"].astype(str)
    else:
        description = pd.Series(df.index.astype(str), index=df.index)

    if "Gene Symbol" in df.columns:
        symbol = df["Gene Symbol"].astype(str)
    elif "Gene.Symbol" in df.columns:
        symbol = df["Gene.Symbol"].astype(str)
    else:
        parsed = description.str.extract(r" GN=(\S+)")[0]
        symbol = parsed.where(parsed.notna() & (parsed.str.strip() != ""), description)

    fullname = _make_unique(description.str.split(" OS=").str[0])

    mode = (label_mode or "fullname").lower()
    if mode == "accession":
        labels = accession
    elif mode == "fullname":
        labels = fullname
    elif mode == "none":
        return pd.Series("", index=df.index)
    else:
        labels = symbol

    blank = labels.isna() | (labels.astype(str).str.strip() == "")
    labels = labels.where(~blank, accession)
    return _make_unique(labels.astype(str))


def build_volcano_figure(
    df: pd.DataFrame,
    pvalue_threshold: float = 0.05,
    fc_threshold: float = 0,
    top_n: int = 10,
    label_mode: str = "fullname",
    test_group_name: str = "Test",
    control_group_name: str = "Control",
) -> go.Figure:
    """
    df: a comparison's full results_df (all kept proteins, not just
        significant ones -- must have logFC and P.Value columns).
    """
    if "logFC" not in df.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="Volcano plot for ANOVA comparisons isn't supported yet.",
            showarrow=False,
            font=dict(size=14),
        )
        fig.update_layout(xaxis_visible=False, yaxis_visible=False)
        return fig

    data = pd.DataFrame(index=df.index)
    data["logFC"] = pd.to_numeric(df["logFC"], errors="coerce")
    data["pvalue"] = pd.to_numeric(df["P.Value"], errors="coerce")
    data["ProteinNames"] = _pick_labels(df, label_mode)
    data = data.dropna(subset=["logFC", "pvalue"])

    up = (data["logFC"] > fc_threshold) & (data["pvalue"] <= pvalue_threshold)
    down = (data["logFC"] < -fc_threshold) & (data["pvalue"] <= pvalue_threshold)
    data["Expression"] = np.select(
        [up, down], ["Up-Regulated", "Down-Regulated"], default="Unchanged"
    )

    fig = go.Figure()
    for expr, color in [
        ("Unchanged", UNCHANGED_COLOR),
        ("Down-Regulated", DOWN_COLOR),
        ("Up-Regulated", UP_COLOR),
    ]:
        subset = data[data["Expression"] == expr]
        if subset.empty:
            continue
        fig.add_trace(
            go.Scattergl(
                x=subset["logFC"],
                y=-np.log10(subset["pvalue"]),
                mode="markers",
                name=expr,
                marker=dict(color=color, size=6, opacity=0.8),
                text=subset["ProteinNames"],
                hovertemplate="%{text}<br>log2FC=%{x:.2f}<br>-log10(p)=%{y:.2f}<extra></extra>",
            )
        )

    if top_n and top_n > 0 and label_mode != "none":
        up_top = (
            data[data["Expression"] == "Up-Regulated"]
            .sort_values("logFC", ascending=False)
            .head(int(top_n))
        )
        down_top = (
            data[data["Expression"] == "Down-Regulated"].sort_values("logFC").head(int(top_n))
        )
        for _, row in pd.concat([up_top, down_top]).iterrows():
            fig.add_annotation(
                x=row["logFC"],
                y=-np.log10(row["pvalue"]),
                text=row["ProteinNames"],
                showarrow=True,
                arrowhead=0,
                ax=0,
                ay=-15,
                font=dict(size=9),
                bgcolor="rgba(255,255,255,0.7)",
            )

    fig.add_hline(y=-np.log10(pvalue_threshold), line_dash="dash", line_color="black")
    if fc_threshold > 0:
        fig.add_vline(x=fc_threshold, line_dash="dash", line_color="black")
        fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="black")

    fig.update_layout(
        title=f"Volcano Plot: {test_group_name} vs. {control_group_name}",
        xaxis_title="Log2 fold change",
        yaxis_title="-Log10 p-value",
        legend_title="",
        template="plotly_white",
        margin=dict(t=50, b=40),
    )
    return fig
