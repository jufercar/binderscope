from __future__ import annotations

import pytest
from Bio.PDB import PDBParser

from binderscope.config import Mutation, Segment
from binderscope.mapping import SegmentMap
from binderscope.reversion import ReversionError, revert_design

MUTATION = Mutation(position=5, wild_type="ALA", reference_position=105)


def segment_map() -> SegmentMap:
    return SegmentMap([Segment(input_start=1, input_end=10, reference_start=101)])


def run(tmp_path, design, native, reference, strategy="truncate"):
    output = tmp_path / "reverted.pdb"
    result = revert_design(
        design_pdb=design,
        output_pdb=output,
        native_target_pdb=native,
        reference_pdb=reference,
        segment_map=segment_map(),
        mutations=[MUTATION],
        native_chain="D",
        reference_chain="D",
        design_target_chain="A",
        design_binder_chain="B",
        strategy=strategy,
    )
    model = PDBParser(QUIET=True).get_structure("out", str(output))[0]
    return result, model


def test_reverted_residue_is_relabelled(tmp_path, design_pdb, native_target_pdb, reference_pdb):
    result, model = run(tmp_path, design_pdb, native_target_pdb, reference_pdb)
    assert result.reverted == [5]
    assert model["A"][5].get_resname() == "ALA"


def test_truncation_leaves_only_rebuildable_atoms(
    tmp_path, design_pdb, native_target_pdb, reference_pdb
):
    # The engineered lysine side chain must go: Rosetta rebuilds the wild-type
    # side chain from backbone + CB with ideal geometry.
    _, model = run(tmp_path, design_pdb, native_target_pdb, reference_pdb)
    names = {atom.get_name() for atom in model["A"][5]}
    assert names <= {"N", "CA", "C", "O", "OXT", "CB"}
    assert not names & {"CG", "CD", "CE", "NZ"}


def test_donor_strategy_transplants_the_reference_rotamer(
    tmp_path, design_pdb, native_target_pdb, reference_pdb
):
    _, model = run(tmp_path, design_pdb, native_target_pdb, reference_pdb, strategy="donor")
    residue = model["A"][5]
    assert residue.get_resname() == "ALA"
    assert "CB" in residue
    assert not {atom.get_name() for atom in residue} & {"CG", "NZ"}


def test_donor_survives_a_reference_in_another_frame(
    tmp_path, design_pdb, native_target_pdb, reference_pdb, reference_pdb_other_frame
):
    # Regression: the reference carries its own coordinate frame, so it needs
    # its own superposition. Reusing the native-to-design transform grafted
    # side chains tens of ångström from their own backbone, without any error.
    def bond_length(reference):
        _, model = run(tmp_path, design_pdb, native_target_pdb, reference, strategy="donor")
        residue = model["A"][5]
        ca = residue["CA"].get_vector().get_array()
        cb = residue["CB"].get_vector().get_array()
        return float(((cb - ca) ** 2).sum() ** 0.5)

    same_frame = bond_length(reference_pdb)
    other_frame = bond_length(reference_pdb_other_frame)
    assert other_frame == pytest.approx(same_frame, abs=1e-3)


def test_donor_preserves_the_reference_rotamer_geometry(
    tmp_path, design_pdb, native_target_pdb, reference_pdb_other_frame
):
    # The point of grafting is to keep the reference's own side-chain geometry.
    # A global-only superposition left the CA-CB bond stretched by up to 0.4 Å,
    # which chi-only relax cannot repair because it moves torsions, not bonds.
    reference = PDBParser(QUIET=True).get_structure("r", str(reference_pdb_other_frame))[0]
    donor = reference["D"][105]
    expected = donor["CB"] - donor["CA"]

    _, model = run(
        tmp_path, design_pdb, native_target_pdb, reference_pdb_other_frame, strategy="donor"
    )
    grafted = model["A"][5]
    assert (grafted["CB"] - grafted["CA"]) == pytest.approx(expected, abs=1e-4)


def test_donor_needs_backbone_anchors(
    tmp_path, design_pdb, native_target_pdb, reference_pdb, monkeypatch
):
    import binderscope.reversion as reversion

    original = reversion._residue_index

    def strip_anchor(chain):
        index = original(chain)
        # Drop N from the donor residue, leaving too few anchors to graft on.
        if 105 in index and "N" in index[105]:
            index[105].detach_child("N")
        return index

    monkeypatch.setattr(reversion, "_residue_index", strip_anchor)
    with pytest.raises(ReversionError, match="without backbone atoms"):
        run(tmp_path, design_pdb, native_target_pdb, reference_pdb, strategy="donor")


def test_donor_reports_the_reference_superposition(
    tmp_path, design_pdb, native_target_pdb, reference_pdb_other_frame
):
    result, _ = run(
        tmp_path, design_pdb, native_target_pdb, reference_pdb_other_frame, strategy="donor"
    )
    assert result.rmsd_reference is not None
    assert result.rmsd_reference == pytest.approx(0.0, abs=1e-2)


def test_truncate_needs_no_reference_superposition(
    tmp_path, design_pdb, native_target_pdb, reference_pdb
):
    result, _ = run(tmp_path, design_pdb, native_target_pdb, reference_pdb)
    assert result.rmsd_reference is None


def test_unmutated_residues_are_untouched(
    tmp_path, design_pdb, native_target_pdb, reference_pdb
):
    _, model = run(tmp_path, design_pdb, native_target_pdb, reference_pdb)
    assert [r.get_resname() for r in model["A"]].count("GLY") == 9


def test_binder_chain_survives_intact(tmp_path, design_pdb, native_target_pdb, reference_pdb):
    result, model = run(tmp_path, design_pdb, native_target_pdb, reference_pdb)
    assert result.binder_length == 5
    assert len(list(model["B"])) == 5
    assert all(r.get_resname() == "ALA" for r in model["B"])


def test_superposition_onto_identical_geometry_is_exact(
    tmp_path, design_pdb, native_target_pdb, reference_pdb
):
    result, _ = run(tmp_path, design_pdb, native_target_pdb, reference_pdb)
    assert result.rmsd_ca == pytest.approx(0.0, abs=1e-6)
    assert result.anchors == 10


def test_hydrogens_and_heteroatoms_are_dropped(
    tmp_path, design_pdb, native_target_pdb, reference_pdb
):
    _, model = run(tmp_path, design_pdb, native_target_pdb, reference_pdb)
    assert all(
        atom.element not in ("H", "D")
        for chain in model
        for residue in chain
        for atom in residue
    )


def test_missing_chain_is_reported_clearly(tmp_path, design_pdb, native_target_pdb, reference_pdb):
    with pytest.raises(ReversionError, match="chain 'Z' not found"):
        revert_design(
            design_pdb=design_pdb,
            output_pdb=tmp_path / "out.pdb",
            native_target_pdb=native_target_pdb,
            reference_pdb=reference_pdb,
            segment_map=segment_map(),
            mutations=[MUTATION],
            native_chain="D",
            reference_chain="D",
            design_target_chain="Z",
            design_binder_chain="B",
        )


def test_too_few_anchors_is_reported_clearly(
    tmp_path, design_pdb, native_target_pdb, reference_pdb
):
    # A map covering positions the design does not have leaves nothing to align on.
    sparse = SegmentMap([Segment(input_start=900, input_end=902, reference_start=1)])
    with pytest.raises(ReversionError, match="Cα anchors"):
        revert_design(
            design_pdb=design_pdb,
            output_pdb=tmp_path / "out.pdb",
            native_target_pdb=native_target_pdb,
            reference_pdb=reference_pdb,
            segment_map=sparse,
            mutations=[MUTATION],
            native_chain="D",
            reference_chain="D",
        )


def test_unknown_strategy_is_rejected(tmp_path, design_pdb, native_target_pdb, reference_pdb):
    with pytest.raises(ValueError, match="unknown reversion strategy"):
        run(tmp_path, design_pdb, native_target_pdb, reference_pdb, strategy="magic")
