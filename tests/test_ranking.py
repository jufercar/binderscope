from __future__ import annotations

import pandas as pd
import pytest

from binderscope.config import Filter, RankMetric
from binderscope.ranking import apply_filters, normalise, rank_designs


def test_normalise_maps_to_unit_range():
    result = normalise(pd.Series([0.0, 5.0, 10.0]), higher_is_better=True)
    assert list(result) == [0.0, 0.5, 1.0]


def test_normalise_inverts_when_lower_is_better():
    result = normalise(pd.Series([0.0, 5.0, 10.0]), higher_is_better=False)
    assert list(result) == [1.0, 0.5, 0.0]


def test_constant_column_carries_no_signal():
    result = normalise(pd.Series([3.0, 3.0, 3.0]))
    assert list(result) == [0.5, 0.5, 0.5]


def test_non_numeric_values_do_not_raise():
    result = normalise(pd.Series(["a", "b"]))
    assert list(result) == [0.5, 0.5]


def test_filters_combine_conjunctively():
    df = pd.DataFrame({"dG": [-30, -5, -40], "sc": [0.8, 0.9, 0.4]})
    filters = [
        Filter(col="dG", label="dG", op="<", value=-10),
        Filter(col="sc", label="sc", op=">", value=0.6),
    ]
    assert list(apply_filters(df, filters)) == [True, False, False]


def test_disabled_and_absent_filters_are_ignored():
    df = pd.DataFrame({"dG": [-30, -5]})
    filters = [
        Filter(col="dG", label="dG", op="<", value=-10, enabled=False),
        Filter(col="missing", label="nope", op=">", value=1),
    ]
    assert list(apply_filters(df, filters)) == [True, True]


METRICS = [
    RankMetric(col="dG", label="dG", higher_is_better=False, weight=2.0),
    RankMetric(col="sc", label="SC", higher_is_better=True, weight=1.0),
]


def test_best_design_ranks_first():
    df = pd.DataFrame(
        {
            "design": ["poor", "best", "middling"],
            "dG": [-10.0, -40.0, -25.0],
            "sc": [0.50, 0.80, 0.65],
        }
    )
    ranked, contributions = rank_designs(df, METRICS)
    assert list(ranked["design"]) == ["best", "middling", "poor"]
    assert list(ranked["rank"]) == [1, 2, 3]
    assert ranked.iloc[0]["score"] == pytest.approx(3.0)
    assert ranked["score_max"].iloc[0] == 3.0
    assert list(contributions.index) == list(ranked.index)


def test_weights_change_the_order():
    df = pd.DataFrame(
        {
            "design": ["energy", "shape"],
            "dG": [-40.0, -10.0],
            "sc": [0.50, 0.90],
        }
    )
    energy_first, _ = rank_designs(
        df,
        [
            RankMetric(col="dG", label="dG", higher_is_better=False, weight=3.0),
            RankMetric(col="sc", label="SC", higher_is_better=True, weight=1.0),
        ],
    )
    shape_first, _ = rank_designs(
        df,
        [
            RankMetric(col="dG", label="dG", higher_is_better=False, weight=1.0),
            RankMetric(col="sc", label="SC", higher_is_better=True, weight=3.0),
        ],
    )
    assert energy_first.iloc[0]["design"] == "energy"
    assert shape_first.iloc[0]["design"] == "shape"


def test_exclusions_apply_before_normalisation():
    # The artefact at dG = +5 would otherwise stretch the scale and compress
    # the two real designs together.
    df = pd.DataFrame(
        {
            "design": ["good", "bad", "artefact"],
            "dG": [-40.0, -20.0, 5.0],
            "sc": [0.80, 0.60, 0.99],
        }
    )
    ranked, _ = rank_designs(
        df, METRICS, exclude=[Filter(col="dG", label="dG", op="<", value=0)]
    )
    assert "artefact" not in list(ranked["design"])
    assert len(ranked) == 2
    assert ranked.iloc[0]["design"] == "good"
    assert ranked.iloc[0]["score"] == pytest.approx(3.0)


def test_missing_metrics_are_reported_but_not_fatal():
    df = pd.DataFrame({"design": ["a", "b"], "dG": [-30.0, -10.0]})
    ranked, contributions = rank_designs(df, METRICS)
    assert ranked.attrs["missing_metrics"] == ["sc"]
    assert list(contributions.columns) == ["dG"]


def test_no_usable_metrics_raises():
    df = pd.DataFrame({"design": ["a"], "other": [1.0]})
    with pytest.raises(ValueError, match="none of the configured ranking metrics"):
        rank_designs(df, METRICS)


def test_empty_metric_list_raises():
    with pytest.raises(ValueError, match="no ranking metrics configured"):
        rank_designs(pd.DataFrame({"dG": [-1.0]}), [])


def test_excluding_everything_yields_an_empty_ranking():
    df = pd.DataFrame({"design": ["a"], "dG": [5.0], "sc": [0.5]})
    ranked, _ = rank_designs(
        df, METRICS, exclude=[Filter(col="dG", label="dG", op="<", value=0)]
    )
    assert ranked.empty
