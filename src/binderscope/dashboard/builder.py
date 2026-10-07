"""Render a self-contained interactive dashboard from a metrics table.

The output is one HTML file with the data inlined, so it can be emailed, put in
a supplement or opened from a USB stick with no server and no install. Plotly is
the only external asset.

Everything the dashboard says about the target is generated from the
configuration, never hard-coded: a run with no ``reference`` section simply has
fewer context panels.
"""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import pandas as pd

from ..config import Config
from ..mapping import SegmentMap

TEMPLATE = Path(__file__).with_name("template.html")

#: Preferred left-to-right column order. Unlisted columns keep their own order
#: at the end, so new metrics show up without touching this list.
COLUMN_ORDER = [
    "design", "model", "binder_len",
    "dG", "dSASA", "dG_dSASA_ratio", "sc", "packstat",
    "n_hbonds", "hbond_pct", "unsat_hbonds", "unsat_hbonds_pct",
    "n_interface_res", "interface_hydro", "interface_sasa_pct",
    "binder_score", "binder_sasa", "surface_hydro",
    "binder_helix_pct", "binder_sheet_pct", "binder_loop_pct",
    "iface_helix_pct", "iface_sheet_pct", "iface_loop_pct",
    "rmsd_reference", "rmsd_design_native", "n_superposition_anchors",
    "clashes_pre_relax",
]

DEFAULT_VISIBLE = [
    "design", "binder_len", "dG", "dSASA", "sc", "packstat",
    "n_hbonds", "unsat_hbonds", "Average_i_pTM", "Average_ipSAE",
    "Average_pLDDT", "rmsd_reference",
]

#: Radar axis labels and directions for metrics this package computes itself.
#: Columns coming from the predictor are added automatically.
AXIS_HINTS: dict[str, tuple[str, bool]] = {
    "dG": ("dG", False),
    "dSASA": ("dSASA", True),
    "sc": ("SC", True),
    "packstat": ("packstat", True),
    "n_hbonds": ("H-bonds", True),
    "unsat_hbonds": ("BUNS", False),
    "interface_sasa_pct": ("iface_SASA%", True),
    "surface_hydro": ("surf_hydro", False),
    "interface_hydro": ("iface_hydro", False),
    "binder_helix_pct": ("helix%", True),
    "binder_score": ("binder_E", False),
    "clashes_pre_relax": ("clashes_pre", False),
    "rmsd_reference": ("RMSD_ref", False),
    "Average_i_pTM": ("i_pTM", True),
    "Average_ipSAE": ("ipSAE", True),
    "Average_pLDDT": ("pLDDT", True),
    "Average_Binder_pLDDT": ("binder_pLDDT", True),
}


def _slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "designs"


def _axis_meta(df: pd.DataFrame, config: Config) -> list[dict]:
    """Radar axes: configured ranking metrics first, then any other known metric."""
    seen: set[str] = set()
    axes: list[dict] = []

    ordered: list[tuple[str, str, bool]] = [
        (m.col, m.label, m.higher_is_better) for m in config.rank_metrics
    ]
    ordered += [
        (col, label, higher)
        for col, (label, higher) in AXIS_HINTS.items()
        if col not in {c for c, _, _ in ordered}
    ]

    for col, label, higher in ordered:
        if col in seen or col not in df.columns:
            continue
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        if values.empty or values.min() == values.max():
            continue
        seen.add(col)
        axes.append(
            {
                "col": col,
                "label": label,
                "higher": bool(higher),
                "min": round(float(values.min()), 4),
                "max": round(float(values.max()), 4),
            }
        )
    return axes


def _table_columns(df: pd.DataFrame) -> list[str]:
    ordered = [c for c in COLUMN_ORDER if c in df.columns]
    spatial = sorted(c for c in df.columns if c.startswith(("clash_", "V_")) and c not in ordered)
    predictor = [
        c for c in df.columns
        if c.startswith(("Average_", "MPNN_")) or c in ("Rank", "Length", "Seed", "Helicity")
    ]
    tail = [
        c for c in df.columns
        if c not in ordered + spatial + predictor and c != "design_base"
    ]
    return ordered + spatial + predictor + tail


def _context_html(config: Config, n_designs: int) -> str:
    """Build the context panels from configuration alone."""
    cards: list[str] = []
    esc = html.escape

    overview = [
        f"<p><b>Designs analysed:</b> {n_designs}</p>",
        f"<p><b>Target chain:</b> {esc(config.target_chain)} &nbsp;·&nbsp; "
        f"<b>Binder chain:</b> {esc(config.binder_chain)}</p>",
    ]
    if config.description:
        overview.insert(0, f"<p>{esc(config.description)}</p>")
    cards.append(_card(f"Run: {esc(config.name)}", "".join(overview)))

    reference = config.reference
    if reference and reference.segments:
        segment_map = SegmentMap.from_reference(reference)
        rows = "".join(
            "<tr>"
            f"<td>{esc(str(row['label']))}</td>"
            f"<td>{esc(str(row['design_range']))}</td>"
            f"<td>{esc(str(row['reference_range']))}</td>"
            f"<td>{esc(str(row['offset']))}</td>"
            f"<td>{row['length']} res</td>"
            "</tr>"
            for row in segment_map.describe()
        )
        linkers = sorted(segment_map.linker_positions)
        body = (
            "<p>The target chain is renumbered from 1 in the design PDBs. "
            "Analysis maps those positions back onto the reference structure "
            f"<code>{esc(reference.structure.name)}</code> "
            f"(chain {esc(reference.target_chain)}).</p>"
            '<table class="ctx-table"><thead><tr>'
            "<th>Fragment</th><th>Design residues</th><th>Reference residues</th>"
            "<th>Offset</th><th>Length</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>"
        )
        if linkers:
            body += (
                f"<p><b>{len(linkers)}</b> linker positions have no reference "
                "equivalent and are excluded from every superposition.</p>"
            )
        body += (
            f"<p><b>{len(segment_map)}</b> mapped residues are available as "
            "superposition anchors.</p>"
        )
        cards.append(_card("Target mapping", body))

    if reference and reference.revert_mutations:
        rows = "".join(
            "<tr>"
            f"<td>{m.position}</td>"
            f"<td>{m.reference_position}</td>"
            f"<td>{esc(m.wild_type)}</td>"
            "</tr>"
            for m in reference.revert_mutations
        )
        cards.append(
            _card(
                "Reverted substitutions",
                "<p>The design target carries engineered substitutions. Each is "
                "restored to wild type before any metric is computed, so the "
                "scored interface is the one that matters biologically.</p>"
                '<table class="ctx-table"><thead><tr><th>Design pos.</th>'
                "<th>Reference pos.</th><th>Wild type</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>",
            )
        )

    steps = [
        "Load the BindCraft complex (target chain + designed binder).",
    ]
    if config.needs_reversion:
        steps += [
            "Superpose the as-designed target onto the design using every mapped "
            "non-linker C&alpha;, then restore wild-type residues in place.",
            "Write the reverted complex: wild-type target in the predicted frame, "
            "binder untouched.",
        ]
    steps += [
        "<b>PyRosetta:</b> chi-only FastRelax &rarr; InterfaceAnalyzerMover &rarr; "
        "BuriedUnsatHbonds (DAlphaBall) &rarr; total energy, SASA and layer metrics.",
    ]
    if reference and reference.clash_targets:
        labels = ", ".join(
            f"<code>{esc(t.label)}</code> (expected {esc(t.expect)})"
            for t in reference.clash_targets
        )
        steps.append(
            "<b>Placement geometry:</b> superpose the complex into the reference "
            f"frame and count binder heavy atoms contacting {labels}."
        )
    if config.bindcraft_stats:
        steps.append(
            "Join the predictor statistics BindCraft already wrote "
            "(<code>final_design_stats.csv</code>)."
        )
    cards.append(
        _card("Pipeline", "<ol>" + "".join(f"<li>{s}</li>" for s in steps) + "</ol>")
    )

    if config.rank_metrics:
        rows = "".join(
            "<tr>"
            f"<td><b>{esc(m.label)}</b></td>"
            f"<td>&times;{m.weight:g}</td>"
            f"<td>{'&uarr;' if m.higher_is_better else '&darr;'}</td>"
            f"<td>{esc(m.description)}</td>"
            "</tr>"
            for m in config.rank_metrics
        )
        note = (
            "<p>Each metric is min-max normalised to [0&ndash;1] in its own good "
            "direction and multiplied by its weight; the score is the sum. "
            f"Maximum possible score: <b>{config.total_weight:g}</b>. "
            "Weights encode priorities and live in the config file.</p>"
        )
        if config.rank_exclude:
            excluded = ", ".join(
                f"<code>{esc(f.col)} {esc(f.op)} {f.value:g}</code>"
                for f in config.rank_exclude
            )
            note += (
                f"<p>Rows are kept for ranking only when {excluded}, so artefacts "
                "do not stretch the normalisation scale.</p>"
            )
        cards.append(
            _card(
                "Multi-criteria ranking",
                note
                + '<table class="ctx-table"><thead><tr><th>Metric</th><th>Weight</th>'
                "<th>Dir.</th><th>Meaning</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>",
            )
        )
    return "\n".join(cards)


def _card(title: str, body: str) -> str:
    return f'  <div class="ctx-card">\n    <h3>{title}</h3>\n    {body}\n  </div>'


def build_dashboard(
    df: pd.DataFrame,
    config: Config,
    output: str | Path | None = None,
) -> Path:
    """Write the dashboard and return its path."""
    if df.empty:
        raise ValueError("cannot build a dashboard from an empty metrics table")

    output = Path(output) if output else config.dashboard_html
    output.parent.mkdir(parents=True, exist_ok=True)

    table_columns = _table_columns(df)
    visible = [c for c in DEFAULT_VISIBLE if c in table_columns]
    visible += [
        c for c in (m.col for m in config.rank_metrics)
        if c in table_columns and c not in visible
    ]
    if not visible:
        visible = table_columns[: min(12, len(table_columns))]

    records = df[table_columns].fillna("").round(4).to_dict(orient="records")
    rank_metrics = [
        {
            "col": m.col,
            "label": m.label,
            "higher": m.higher_is_better,
            "w": m.weight,
        }
        for m in config.rank_metrics
        if m.col in df.columns
    ]
    filters = [
        {
            "col": f.col,
            "label": f.label,
            "op": f.op,
            "val": f.value,
            "enabled": f.enabled,
        }
        for f in config.filters
        if f.col in df.columns
    ]
    exclude = [
        {"col": f.col, "op": f.op, "val": f.value}
        for f in config.rank_exclude
        if f.col in df.columns
    ]

    subtitle = (
        f"{len(records)} designs · run {config.name} · "
        "PyRosetta interface metrics"
        + (" + placement geometry" if config.reference else "")
        + (" + predictor statistics" if config.bindcraft_stats else "")
    )

    html_out = TEMPLATE.read_text()
    for placeholder, value in {
        "__TITLE__": html.escape(f"{config.name} — binder design triage"),
        "__HSTXT__": html.escape(subtitle),
        "__CONTEXT_HTML__": _context_html(config, len(records)),
        "__DATA__": json.dumps(records),
        "__RMETA__": json.dumps(_axis_meta(df, config)),
        "__FDEFS__": json.dumps(filters),
        "__TCOLS__": json.dumps(table_columns),
        "__DVCOLS__": json.dumps(visible),
        "__RANK_METRICS__": json.dumps(rank_metrics),
        "__RANK_EXCLUDE__": json.dumps(exclude),
        "__SLUG__": _slug(config.name),
    }.items():
        html_out = html_out.replace(placeholder, value)

    leftover = re.findall(r"__[A-Z_]+__", html_out)
    if leftover:
        raise RuntimeError(f"unfilled dashboard placeholders: {sorted(set(leftover))}")

    output.write_text(html_out)
    return output
