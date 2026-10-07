from __future__ import annotations

import pytest

from binderscope.config import ConfigError, Segment
from binderscope.mapping import SegmentMap


def make_map() -> SegmentMap:
    """Three segments joined by two linkers, as in a spliced construct."""
    segments = [
        Segment(input_start=1, input_end=73, reference_start=81, label="seg1"),
        Segment(input_start=77, input_end=83, reference_start=417, label="seg2"),
        Segment(input_start=86, input_end=141, reference_start=435, label="seg3"),
    ]
    linkers = set(range(74, 77)) | set(range(84, 86))
    return SegmentMap(segments, linkers)


def test_offsets_are_derived_not_declared():
    segment = Segment(input_start=1, input_end=73, reference_start=81)
    assert segment.offset == 80
    assert segment.length == 73


def test_forward_mapping_across_segments():
    mapping = make_map()
    assert mapping.to_reference(1) == 81
    assert mapping.to_reference(73) == 153
    assert mapping.to_reference(77) == 417
    assert mapping.to_reference(83) == 423
    assert mapping.to_reference(86) == 435
    assert mapping.to_reference(141) == 490


def test_linkers_and_out_of_range_have_no_reference():
    mapping = make_map()
    for position in (74, 75, 76, 84, 85):
        assert mapping.to_reference(position) is None
        assert mapping.is_linker(position)
    assert mapping.to_reference(0) is None
    assert mapping.to_reference(999) is None


def test_reverse_mapping_round_trips():
    mapping = make_map()
    for design_pos, reference_pos in mapping.pairs:
        assert mapping.to_design(reference_pos) == design_pos


def test_anchor_count_matches_non_linker_residues():
    mapping = make_map()
    # 73 + 7 + 56 residues, with the linkers excluded by construction
    assert len(mapping) == 136
    assert len(mapping.design_positions) == 136


def test_overlapping_segments_are_rejected():
    segments = [
        Segment(input_start=1, input_end=50, reference_start=1),
        Segment(input_start=40, input_end=60, reference_start=200),
    ]
    with pytest.raises(ValueError, match="more than one segment"):
        SegmentMap(segments)


def test_reversed_range_is_rejected():
    with pytest.raises(ConfigError, match="input_end < input_start"):
        Segment(input_start=50, input_end=10, reference_start=1)


def test_describe_reports_derived_ranges():
    rows = make_map().describe()
    assert rows[0]["reference_range"] == "81–153"
    assert rows[0]["offset"] == "+80"
    assert rows[1]["length"] == 7
