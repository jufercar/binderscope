"""Placement geometry in a reference frame.

These exercise the module that answers "did the binder end up where it was
meant to?", which no interface energy can tell you.
"""

from __future__ import annotations

import pytest
from conftest import BINDER, TARGET_ENGINEERED, write_pdb

from binderscope.config import ClashTarget, GeometryOptions, Reference, Segment
from binderscope.mapping import SegmentMap
from binderscope.spatial import ReferenceFrame, SpatialError


def make_reference(structure, *, volume=False, chains=("D", "B")) -> Reference:
    target, partner = chains
    return Reference(
        structure=structure,
        target_chain=target,
        segments=[Segment(input_start=1, input_end=10, reference_start=101)],
        clash_targets=[
            ClashTarget(chain=target, label="target", expect="low"),
            ClashTarget(
                chain=partner, label="native_partner", expect="high", measure_volume=volume
            ),
        ],
    )


@pytest.fixture
def frame(reference_pdb):
    reference = make_reference(reference_pdb, volume=True)
    return ReferenceFrame(reference, SegmentMap.from_reference(reference), GeometryOptions())


def test_column_names_follow_the_configured_labels(frame, design_pdb):
    # Downstream configs reference these names, so they are part of the contract.
    result = frame.measure(design_pdb, "A", "B")
    assert set(result) == {
        "rmsd_reference",
        "n_superposition_anchors",
        "clash_target_hard",
        "clash_target_soft",
        "clash_native_partner_hard",
        "clash_native_partner_soft",
        "V_binder_A3",
        "V_inter_native_partner_A3",
        "V_inter_native_partner_pct",
    }


def test_binder_volume_is_reported_once_per_run(reference_pdb, design_pdb):
    # It describes the binder, not the chain it is compared against, so asking
    # for two overlap volumes must not produce two differing binder volumes.
    reference = Reference(
        structure=reference_pdb,
        target_chain="D",
        segments=[Segment(input_start=1, input_end=10, reference_start=101)],
        clash_targets=[
            ClashTarget(chain="D", label="target", measure_volume=True),
            ClashTarget(chain="B", label="native_partner", measure_volume=True),
        ],
    )
    frame = ReferenceFrame(reference, SegmentMap.from_reference(reference))
    result = frame.measure(design_pdb, "A", "B")
    assert "V_inter_target_A3" in result
    assert "V_inter_native_partner_A3" in result
    assert len([k for k in result if k.startswith("V_binder")]) == 1


def test_superposition_uses_every_mapped_residue(frame, design_pdb):
    result = frame.measure(design_pdb, "A", "B")
    assert result["n_superposition_anchors"] == 10


def test_identical_geometry_superposes_exactly(frame, design_pdb):
    assert frame.measure(design_pdb, "A", "B")["rmsd_reference"] == pytest.approx(0.0, abs=1e-6)


def test_a_reference_in_another_frame_still_superposes(
    reference_pdb_other_frame, design_pdb
):
    # The whole point of superposing is that the reference's own frame is
    # irrelevant; a rotated, translated reference must give the same answer.
    reference = make_reference(reference_pdb_other_frame, volume=True)
    frame = ReferenceFrame(reference, SegmentMap.from_reference(reference))
    result = frame.measure(design_pdb, "A", "B")
    assert result["rmsd_reference"] == pytest.approx(0.0, abs=1e-2)


def test_clash_counts_respond_to_the_configured_thresholds(reference_pdb, design_pdb):
    reference = make_reference(reference_pdb)
    mapping = SegmentMap.from_reference(reference)
    tight = ReferenceFrame(reference, mapping, GeometryOptions(clash_hard=0.1, clash_soft=0.2))
    loose = ReferenceFrame(reference, mapping, GeometryOptions(clash_hard=20.0, clash_soft=25.0))
    # Measured against the target chain, which the fixture places 12 Å from the
    # binder, so the count genuinely moves with the threshold.
    assert tight.measure(design_pdb, "A", "B")["clash_target_hard"] == 0
    assert loose.measure(design_pdb, "A", "B")["clash_target_hard"] > 0


def test_clash_count_cannot_exceed_the_binder_atom_count(frame, design_pdb):
    # Atoms are counted, not pairs, which is what keeps the number comparable
    # between designs of different size.
    from Bio.PDB import PDBParser

    from binderscope.geometry import heavy_atoms

    model = PDBParser(QUIET=True).get_structure("d", str(design_pdb))[0]
    n_atoms = len(heavy_atoms(model["B"]))
    result = frame.measure(design_pdb, "A", "B")
    assert result["clash_native_partner_soft"] <= n_atoms


def test_missing_reference_target_chain_is_reported(reference_pdb):
    reference = make_reference(reference_pdb, chains=("Q", "B"))
    with pytest.raises(SpatialError, match="reference chain 'Q' not found"):
        ReferenceFrame(reference, SegmentMap.from_reference(reference))


def test_missing_clash_target_chain_is_reported(reference_pdb):
    reference = Reference(
        structure=reference_pdb,
        target_chain="D",
        segments=[Segment(input_start=1, input_end=10, reference_start=101)],
        clash_targets=[ClashTarget(chain="Z", label="ghost")],
    )
    with pytest.raises(SpatialError, match="clash target chain 'Z' not found"):
        ReferenceFrame(reference, SegmentMap.from_reference(reference))


def test_missing_design_chain_is_reported(frame, design_pdb):
    with pytest.raises(SpatialError, match="chain 'Z' not found"):
        frame.measure(design_pdb, "Z", "B")


def test_unmappable_design_is_reported(reference_pdb, tmp_path):
    # A design whose numbering the segment map does not cover leaves nothing to
    # align on; better to say so than to return a meaningless RMSD.
    reference = make_reference(reference_pdb)
    mapping = SegmentMap([Segment(input_start=900, input_end=910, reference_start=101)])
    frame = ReferenceFrame(reference, mapping)
    design = write_pdb(tmp_path / "odd.pdb", {"A": TARGET_ENGINEERED, "B": BINDER})
    with pytest.raises(SpatialError, match="Cα anchors"):
        frame.measure(design, "A", "B")


def test_a_distant_binder_clashes_with_nothing(reference_pdb, tmp_path):
    reference = make_reference(reference_pdb, volume=True)
    frame = ReferenceFrame(reference, SegmentMap.from_reference(reference))
    far_binder = [(n, "ALA") for n in range(1, 6)]
    design = write_pdb(tmp_path / "far.pdb", {"A": TARGET_ENGINEERED, "B": far_binder})

    # Push the binder far away by rewriting its coordinates.
    lines = []
    for line in design.read_text().splitlines(keepends=True):
        if line.startswith("ATOM") and line[21] == "B":
            x = float(line[30:38]) + 500.0
            line = f"{line[:30]}{x:8.3f}{line[38:]}"
        lines.append(line)
    design.write_text("".join(lines))

    result = frame.measure(design, "A", "B")
    assert result["clash_target_hard"] == 0
    assert result["clash_native_partner_hard"] == 0
    assert result["V_inter_native_partner_pct"] == 0.0
    assert result["V_binder_A3"] > 0


def test_reference_is_not_moved_by_measuring(frame, design_pdb, reference_pdb):
    # The reference coordinates are loaded once and reused across designs; if
    # measuring mutated them, results would drift design by design.
    first = frame.measure(design_pdb, "A", "B")
    second = frame.measure(design_pdb, "A", "B")
    assert first == second
