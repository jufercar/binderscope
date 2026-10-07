"""Residue-numbering map between a design-numbered target and a reference structure.

BindCraft renumbers the target chain from 1, and an engineered target may be a
chimera of several reference fragments joined by linkers. This module turns a
list of configured segments into a bidirectional residue map, so no offset is
ever hard-coded in analysis code.
"""

from __future__ import annotations

from collections.abc import Iterator

from .config import Reference, Segment


class SegmentMap:
    """Maps design residue numbers to reference residue numbers.

    Positions inside a linker, or outside every segment, have no reference
    equivalent and map to ``None``.
    """

    def __init__(self, segments: list[Segment], linkers: set[int] | None = None):
        self._segments = list(segments)
        self._linkers = set(linkers or ())
        self._forward: dict[int, int] = {}
        for seg in self._segments:
            for pos in range(seg.input_start, seg.input_end + 1):
                if pos in self._linkers:
                    continue
                if pos in self._forward:
                    raise ValueError(
                        f"design residue {pos} is claimed by more than one segment"
                    )
                self._forward[pos] = pos + seg.offset
        self._reverse = {v: k for k, v in self._forward.items()}

    @classmethod
    def from_reference(cls, reference: Reference) -> SegmentMap:
        return cls(reference.segments, reference.linker_positions)

    # ── lookups ──────────────────────────────────────────────────────────────
    def to_reference(self, design_position: int) -> int | None:
        """Reference residue number for a design position, or None if unmapped."""
        return self._forward.get(design_position)

    def to_design(self, reference_position: int) -> int | None:
        """Design residue number for a reference position, or None if unmapped."""
        return self._reverse.get(reference_position)

    def is_linker(self, design_position: int) -> bool:
        return design_position in self._linkers

    # ── bulk access ──────────────────────────────────────────────────────────
    @property
    def design_positions(self) -> list[int]:
        return sorted(self._forward)

    @property
    def pairs(self) -> list[tuple[int, int]]:
        """(design_position, reference_position) pairs, sorted by design position."""
        return sorted(self._forward.items())

    @property
    def linker_positions(self) -> set[int]:
        return set(self._linkers)

    @property
    def segments(self) -> list[Segment]:
        return list(self._segments)

    def __len__(self) -> int:
        return len(self._forward)

    def __iter__(self) -> Iterator[tuple[int, int]]:
        return iter(self.pairs)

    def __repr__(self) -> str:
        return (
            f"SegmentMap({len(self._segments)} segments, {len(self._forward)} "
            f"mapped residues, {len(self._linkers)} linker positions)"
        )

    def describe(self) -> list[dict[str, object]]:
        """Per-segment summary, used by the dashboard context panel."""
        rows: list[dict[str, object]] = []
        for seg in self._segments:
            rows.append(
                {
                    "label": seg.label,
                    "design_range": f"{seg.input_start}–{seg.input_end}",
                    "reference_range": (
                        f"{seg.reference_start}–{seg.reference_start + seg.length - 1}"
                    ),
                    "offset": f"{seg.offset:+d}",
                    "length": seg.length,
                }
            )
        return rows
