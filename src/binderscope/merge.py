"""Join structural metrics with the predictor statistics BindCraft already wrote."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

#: Columns taken from ``final_design_stats.csv`` when present. Anything absent is
#: skipped silently, so this works across BindCraft versions.
PREDICTOR_COLUMNS = [
    "Design", "Rank", "Length", "Seed", "Helicity",
    "MPNN_score", "MPNN_seq_recovery", "Sequence",
    "Average_pLDDT", "Average_pTM", "Average_i_pTM",
    "Average_pAE", "Average_i_pAE", "Average_ipSAE",
    "Average_i_pLDDT", "Average_ss_pLDDT",
    "Average_Binder_pLDDT", "Average_Binder_pTM", "Average_Binder_pAE",
    "Average_Hotspot_RMSD", "Average_Target_RMSD", "Average_Binder_RMSD",
    "Average_Relaxed_Clashes", "Average_Unrelaxed_Clashes",
    "Average_Binder_Helix%", "Average_Binder_BetaSheet%", "Average_Binder_Loop%",
    "Average_Interface_Helix%", "Average_Interface_BetaSheet%",
    "Average_Interface_Loop%",
    "Average_ShapeComplementarity", "Average_PackStat",
    "Average_dG", "Average_dSASA",
    "Average_Surface_Hydrophobicity", "Average_Binder_Energy_Score",
    "Average_n_InterfaceResidues", "Average_n_InterfaceHbonds",
    "Average_n_InterfaceUnsatHbonds",
]

_MODEL_SUFFIX = re.compile(r"^(?P<base>.*?)(?:_model(?P<model>\d+))?$")


def split_design_name(name: str) -> tuple[str, int | None]:
    """Split ``<design>_model<N>`` into its base name and model number.

    BindCraft writes one PDB per predicted model, while ``final_design_stats.csv``
    holds one row per design, so the suffix has to come off before joining.
    """
    match = _MODEL_SUFFIX.match(name)
    if not match:
        return name, None
    model = match.group("model")
    return match.group("base"), int(model) if model else None


def merge_predictor_stats(
    metrics: pd.DataFrame,
    stats_csv: str | Path | None,
) -> pd.DataFrame:
    """Left-join BindCraft's per-design statistics onto a structural metrics table.

    Returns ``metrics`` unchanged when no statistics file is configured or the
    file holds no rows.
    """
    if stats_csv is None:
        return metrics
    path = Path(stats_csv)
    if not path.is_file():
        raise FileNotFoundError(f"BindCraft statistics file not found: {path}")

    stats = pd.read_csv(path)
    if stats.empty or "Design" not in stats.columns:
        return metrics

    present = [c for c in PREDICTOR_COLUMNS if c in stats.columns]
    subset = stats[present].copy()
    subset["design_base"] = subset["Design"]
    subset = subset.drop(columns=["Design"])

    overlapping = (set(subset.columns) & set(metrics.columns)) - {"design_base"}
    if overlapping:
        subset = subset.drop(columns=sorted(overlapping))

    return metrics.merge(subset, on="design_base", how="left")
