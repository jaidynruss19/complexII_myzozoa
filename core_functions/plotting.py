from pathlib import Path
import pandas as pd
import altair as alt

alt.data_transformers.disable_max_rows()disable_max_rows()


def save_chart(chart: "alt.Chart", out_html: Path) -> Path:
    """Save an Altair chart to HTML and return the output path."""
    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    chart.save(str(out_html))
    return out_html


def make_dotplot(
    plot_data: pd.DataFrame,
    taxname_order: list[str] | None = None,
    y_order: list[str] | None = None,
    class_order: list[str] | None = None,
    title: str = "MMseqs2 hits dotplot",
) -> "alt.Chart":

    required = {"taxname", "subunit", "class_group", "identity_bin"}
    missing = required - set(plot_data.columns)
    if missing:
        raise ValueError(f"plot_data missing required columns: {sorted(missing)}")

    df = plot_data.copy()

    # Default orders (if you already computed these in notebook, pass them in)
    if taxname_order is None:
        taxname_order = df["taxname"].dropna().astype(str).drop_duplicates().tolist()
    if y_order is None:
        y_order = df["subunit"].dropna().astype(str).drop_duplicates().tolist()
    if class_order is None:
        class_order = df["class_group"].dropna().astype(str).drop_duplicates().tolist()

    # Width heuristic taken from your notebook pattern
    fig_width = min(1600, 18 * len(taxname_order))

    # Opacity bins taken from notebook
    opacity_domain = ["0–25% (low)", "25–60% (mid)", "60–100% (high)"]
    opacity_range = [0.25, 0.7, 1.0]

    # Tooltips: keep it robust if some columns don't exist
    tooltip_cols = []
    for col in ["taxid", "taxname", "class_group", "subunit"]:
        if col in df.columns:
            tooltip_cols.append(f"{col}:N")

    # numeric-ish tooltips (only if present)
    if "max_ident" in df.columns:
        tooltip_cols.append(alt.Tooltip("max_ident:Q", format=".1f", title="Max % identity"))
    if "evalue" in df.columns:
        tooltip_cols.append(alt.Tooltip("evalue:Q", format=".2e"))
    if "qcov" in df.columns:
        tooltip_cols.append(alt.Tooltip("qcov:Q", format=".2f"))
    if "tcov" in df.columns:
        tooltip_cols.append(alt.Tooltip("tcov:Q", format=".2f"))

    chart = (
        alt.Chart(df)
        .mark_circle(size=150)
        .encode(
            x=alt.X(
                "taxname:N",
                sort=taxname_order,
                axis=alt.Axis(labelAngle=-90, labelFontSize=14, titleFontSize=14),
                title="Species (tree order)",
            ),
            y=alt.Y(
                "subunit:N",
                sort=y_order,
                axis=alt.Axis(labelFontSize=13, titleFontSize=14),
                title="Complex II subunit",
            ),
            color=alt.Color(
                "class_group:N",
                title="Taxonomic Class",
                scale=alt.Scale(scheme="tableau10"),
                sort=class_order,
            ),
            opacity=alt.Opacity(
                "identity_bin:N",
                scale=alt.Scale(domain=opacity_domain, range=opacity_range),
                title="% Id",
            ),
            tooltip=tooltip_cols if tooltip_cols else None,
        )
        .properties(title=title, width=fig_width)
    )

    return chart


def save_standard_plots(
    plot_data: pd.DataFrame,
    out_dir: Path,
    taxname_order: list[str] | None = None,
    y_order: list[str] | None = None,
    class_order: list[str] | None = None,
) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dot = make_dotplot(
        plot_data,
        taxname_order=taxname_order,
        y_order=y_order,
        class_order=class_order,
        title="MMseqs2 hits dotplot",
    )
    out_dot = save_chart(dot, out_dir / "dotplot_mmseqs_hits.html")

    return {"dotplot": out_dot}
