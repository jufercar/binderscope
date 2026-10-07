from __future__ import annotations

import numpy as np

from binderscope.geometry import count_contacting_atoms, volume_overlap


def test_no_contacts_when_far_apart():
    a = np.zeros((5, 3), dtype=np.float32)
    b = np.full((5, 3), 100.0, dtype=np.float32)
    assert count_contacting_atoms(a, b, 2.0) == 0


def test_contacts_count_query_atoms_not_pairs():
    # One query atom sitting among many neighbours still counts once, which is
    # what makes the metric comparable between designs.
    a = np.zeros((1, 3), dtype=np.float32)
    b = np.array([[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]], dtype=np.float32)
    assert count_contacting_atoms(a, b, 2.0) == 1


def test_contacts_respect_the_threshold():
    a = np.array([[0.0, 0, 0]], dtype=np.float32)
    b = np.array([[2.2, 0, 0]], dtype=np.float32)
    assert count_contacting_atoms(a, b, 2.0) == 0
    assert count_contacting_atoms(a, b, 2.5) == 1


def test_empty_input_is_not_an_error():
    empty = np.zeros((0, 3), dtype=np.float32)
    some = np.zeros((3, 3), dtype=np.float32)
    assert count_contacting_atoms(empty, some, 2.0) == 0
    assert count_contacting_atoms(some, empty, 2.0) == 0


def test_identical_selections_overlap_completely():
    coords = np.array(
        [[0, 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]], dtype=np.float32
    )
    result = volume_overlap(coords, coords)
    assert result["volume"] > 0
    assert result["intersection_pct"] == 100.0


def test_disjoint_selections_do_not_overlap():
    a = np.zeros((4, 3), dtype=np.float32)
    b = np.full((4, 3), 60.0, dtype=np.float32)
    result = volume_overlap(a, b)
    assert result["volume"] > 0
    assert result["intersection"] == 0.0
    assert result["intersection_pct"] == 0.0


def test_partial_overlap_is_between_the_extremes():
    a = np.array([[0, 0, 0], [2, 0, 0], [4, 0, 0]], dtype=np.float32)
    b = np.array([[4, 0, 0], [6, 0, 0], [8, 0, 0]], dtype=np.float32)
    pct = volume_overlap(a, b)["intersection_pct"]
    assert 0.0 < pct < 100.0


def test_empty_query_reports_zero_volume():
    empty = np.zeros((0, 3), dtype=np.float32)
    result = volume_overlap(empty, np.zeros((3, 3), dtype=np.float32))
    assert result == {"volume": 0.0, "intersection": 0.0, "intersection_pct": 0.0}
