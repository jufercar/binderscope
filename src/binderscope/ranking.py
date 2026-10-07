"""Multi-criteria ranking of designs.

No single metric picks a binder. Interface energy, shape complementarity,
predictor confidence and placement geometry each capture something the others
miss, and they trade off against each other. This module normalises every
configured metric to [0, 1] in its own good direction, multiplies by a weight,
and sums. The per-metric contributions are kept alongside the total, because the
reason a design ranks where it does matters more than the rank itself.

Weights are a declaration of priorities, not a discovered truth: they belong in
the config file where a reader can see and change them.
"""

from __future__ import annotations

import pandas as pd

from .config import Filter, RankMetric


def normalise(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    """Min-max scale to [0, 1], flipped when lower values are better.

    A constant column carries no information, so it maps to 0.5 throughout
    rather than to an arbitrary extreme.
    """
    values = pd.to_numeric(series, errors="coerce")
    low, high = values.min(), values.max()
    if pd.isna(low) or pd.isna(high) or high == low:
        return pd.Series(0.5, index=series.index)
    scaled = (values - low) / (high - low)
    return scaled if higher_is_better else 1 - scaled


def apply_filters(df: pd.DataFrame, filters: list[Filter]) -> pd.Series:
    """Boolean mask of rows passing every enabled filter."""
    mask = pd.Series(True, index=df.index)
    for flt in filters:
        if not flt.enabled or flt.col not in df.columns:
            continue
        values = pd.to_numeric(df[flt.col], errors="coerce")
        if flt.op == "<":
            mask &= values < flt.value
        elif flt.op == ">":
            mask &= values > flt.value
        else:
            mask &= values == flt.value
    return mask


def rank_designs(
    df: pd.DataFrame,
    metrics: list[RankMetric],
    exclude: list[Filter] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score and sort designs.

    Returns ``(ranked, contributions)``: the input rows with ``score`` and
    ``rank`` columns added, sorted best first, and the matching per-metric
    weighted contributions.

    ``exclude`` drops rows before normalisation — important, because a single
    artefact (a positive interface energy, say) would otherwise stretch the
    scale and compress every real design into a narrow band.
    """
    if not metrics:
        raise ValueError("no ranking metrics configured")

    working = df.copy()
    if exclude:
        working = working[apply_filters(working, exclude)].copy()
    if working.empty:
        empty = pd.DataFrame(index=working.index)
        return working.assign(score=pd.Series(dtype=float), rank=pd.Series(dtype=int)), empty

    usable = [m for m in metrics if m.col in working.columns]
    missing = [m.col for m in metrics if m.col not in working.columns]
    if not usable:
        raise ValueError(
            "none of the configured ranking metrics are present in the data: "
            + ", ".join(m.col for m in metrics)
        )

    contributions = pd.DataFrame(index=working.index)
    for metric in usable:
        contributions[metric.label] = (
            normalise(working[metric.col], metric.higher_is_better) * metric.weight
        )

    working["score"] = contributions.sum(axis=1)
    working["score_max"] = sum(m.weight for m in usable)
    working = working.sort_values("score", ascending=False)
    working["rank"] = range(1, len(working) + 1)

    if missing:
        working.attrs["missing_metrics"] = missing
    return working, contributions.loc[working.index]
