"""Revert engineered point mutations in a design's target chain to wild type.

Design campaigns often run against a modified target: surface residues swapped
to steer the binder away from an unwanted face, fragments spliced together, and
so on. Scoring the design against that modified target measures the wrong
interface. This module rebuilds the wild-type target in the design's own
coordinate frame, keeping the designed binder untouched.

Two strategies are available:

``truncate`` (default)
    Keep the backbone and ``CB``, relabel the residue as wild type, and let
    Rosetta rebuild the side chain with ideal geometry before the chi-only
    relax repacks it. Chemically consistent for any substitution.

``donor``
    Transplant heavy-atom coordinates from the reference structure. Preserves
    the experimentally observed rotamer, but requires the reference residue to
    be present and complete. The reference is superposed onto the design
    through the residue map before anything is grafted, so it does not matter
    what frame the reference file happens to be in.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

from .config import Mutation
from .mapping import SegmentMap

BACKBONE = ("N", "CA", "C", "O", "OXT")
KEEP_ON_TRUNCATE = BACKBONE + ("CB",)


class ReversionError(RuntimeError):
    """Raised when a design cannot be reverted (missing chains or no anchor atoms)."""


@dataclass
class ReversionResult:
    rmsd_ca: float
    binder_length: int
    anchors: int
    reverted: list[int]
    #: Cα RMSD of the reference onto the design frame. Only the ``donor``
    #: strategy needs that superposition, so it is None otherwise.
    rmsd_reference: float | None = None


def _residue_index(chain) -> dict[int, object]:
    return {r.get_id()[1]: r for r in chain if r.get_id()[0] == " "}


def _superpose(fixed: list, moving: list, what: str, moved_atoms: list) -> float:
    """Superpose ``moving`` onto ``fixed`` and apply the result to ``moved_atoms``."""
    if len(fixed) < 3:
        raise ReversionError(
            f"only {len(fixed)} Cα anchors available to superpose {what} "
            "(need at least 3); check the segment map and chain ids"
        )
    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    superimposer.apply(moved_atoms)
    return float(superimposer.rms)


def revert_design(
    design_pdb: str | Path,
    output_pdb: str | Path,
    native_target_pdb: str | Path,
    reference_pdb: str | Path,
    segment_map: SegmentMap,
    mutations: list[Mutation],
    native_chain: str,
    reference_chain: str,
    design_target_chain: str = "A",
    design_binder_chain: str = "B",
    strategy: str = "truncate",
) -> ReversionResult:
    """Write ``output_pdb`` holding the wild-type target plus the original binder.

    The native target is superposed onto the design's target chain using the
    Cα atoms of every mapped, non-linker position, so the rebuilt target lands
    in the design's (predicted) frame rather than the reference's.
    """
    if strategy not in ("truncate", "donor"):
        raise ValueError(f"unknown reversion strategy: {strategy!r}")

    parser = PDBParser(QUIET=True)
    design = parser.get_structure("design", str(design_pdb))[0]
    native = parser.get_structure("native", str(native_target_pdb))[0]
    reference = parser.get_structure("reference", str(reference_pdb))[0]

    for chain_id, model, what in (
        (design_target_chain, design, "design target"),
        (design_binder_chain, design, "design binder"),
        (native_chain, native, "native target"),
        (reference_chain, reference, "reference"),
    ):
        if chain_id not in model:
            raise ReversionError(f"chain '{chain_id}' not found ({what})")

    native_target = native[native_chain]
    design_target = design[design_target_chain]
    design_binder = design[design_binder_chain]
    reference_target = reference[reference_chain]

    design_index = _residue_index(design_target)
    native_index = _residue_index(native_target)
    reference_index = _residue_index(reference_target)

    # ── superpose native target onto the design's target chain ───────────────
    fixed, moving = [], []
    for position in segment_map.design_positions:
        native_res = native_index.get(position)
        design_res = design_index.get(position)
        if native_res is None or design_res is None:
            continue
        if "CA" in native_res and "CA" in design_res:
            fixed.append(design_res["CA"])
            moving.append(native_res["CA"])

    anchors = len(fixed)
    rmsd_ca = _superpose(
        fixed, moving, "the native target onto the design", list(native_target.get_atoms())
    )

    # ── place the reference in the design frame, for donor transplants ───────
    # The reference is an independent structure: its coordinates are in its own
    # frame, which is NOT the native target's. Reusing the native-to-design
    # transform here would place donor side chains tens of ångström away from
    # the backbone they are grafted onto, silently. So the reference gets its
    # own superposition, through the same residue mapping.
    rmsd_reference: float | None = None
    if strategy == "donor":
        ref_fixed, ref_moving = [], []
        for design_pos, reference_pos in segment_map.pairs:
            design_res = design_index.get(design_pos)
            reference_res = reference_index.get(reference_pos)
            if design_res is None or reference_res is None:
                continue
            if "CA" in design_res and "CA" in reference_res:
                ref_fixed.append(design_res["CA"])
                ref_moving.append(reference_res["CA"])

        rmsd_reference = _superpose(
            ref_fixed,
            ref_moving,
            "the reference onto the design",
            list(reference_target.get_atoms()),
        )

    # ── build the reverted target chain ──────────────────────────────────────
    by_position = {m.position: m for m in mutations}
    out_target = Chain(design_target_chain)
    reverted: list[int] = []

    for position in sorted(native_index):
        residue = copy.deepcopy(native_index[position])
        mutation = by_position.get(position)

        if mutation is not None:
            if strategy == "donor":
                donor = reference_index.get(mutation.reference_position)
                if donor is None:
                    raise ReversionError(
                        f"reference residue {mutation.reference_position} missing; "
                        "cannot transplant side chain (use strategy 'truncate')"
                    )
                residue = _transplant(residue, donor, mutation)
            else:
                residue = _truncate(residue, mutation)
            reverted.append(position)

        residue.detach_parent()
        out_target.add(residue)

    # ── binder chain, untouched ──────────────────────────────────────────────
    out_binder = Chain(design_binder_chain)
    binder_length = 0
    for residue in design_binder:
        if residue.get_id()[0] != " ":
            continue
        clone = copy.deepcopy(residue)
        clone.detach_parent()
        out_binder.add(clone)
        binder_length += 1

    structure = Structure("reverted")
    model = Model(0)
    model.add(out_target)
    model.add(out_binder)
    structure.add(model)

    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), select=_NoHydrogens())

    return ReversionResult(
        rmsd_ca=rmsd_ca,
        binder_length=binder_length,
        anchors=anchors,
        reverted=reverted,
        rmsd_reference=rmsd_reference,
    )


def _truncate(residue, mutation: Mutation):
    """Strip to backbone + CB and relabel; Rosetta rebuilds the side chain."""
    for atom in [a for a in residue if a.get_name() not in KEEP_ON_TRUNCATE]:
        residue.detach_child(atom.get_id())
    if mutation.wild_type == "GLY" and "CB" in residue:
        residue.detach_child("CB")
    residue.resname = mutation.wild_type
    return residue


GRAFT_ANCHORS = ("N", "CA", "C")


def _transplant(residue, donor, mutation: Mutation):
    """Replace side-chain atoms with the reference rotamer.

    The donor is superposed onto this residue by its own backbone N-CA-C before
    anything is copied. A global superposition of the whole reference chain is
    not enough: it is accurate to a few tenths of an ångström overall, but at
    any individual residue the two backbones still differ, and grafting across
    that gap stretches the CA-CB bond. Chi-only relax rebuilds torsions, not
    bond lengths, so such distortion would survive into the scored pose.
    """
    donor = copy.deepcopy(donor)

    missing = [name for name in GRAFT_ANCHORS if name not in residue or name not in donor]
    if missing:
        raise ReversionError(
            f"residue {mutation.position}: cannot graft a side chain without "
            f"backbone atoms {', '.join(missing)} in both the target and the "
            "reference (use strategy 'truncate')"
        )

    _superpose(
        [residue[name] for name in GRAFT_ANCHORS],
        [donor[name] for name in GRAFT_ANCHORS],
        f"the reference rotamer at residue {mutation.position}",
        list(donor.get_atoms()),
    )

    for atom in [a for a in residue if a.get_name() not in BACKBONE]:
        residue.detach_child(atom.get_id())
    for atom in donor:
        name = atom.get_name()
        if name in BACKBONE or atom.element in ("H", "D"):
            continue
        clone = copy.deepcopy(atom)
        clone.detach_parent()
        residue.add(clone)
    residue.resname = mutation.wild_type
    return residue


class _NoHydrogens:
    """PDBIO selector that drops hydrogens and heteroatoms."""

    def accept_model(self, model) -> int:
        return 1

    def accept_chain(self, chain) -> int:
        return 1

    def accept_residue(self, residue) -> int:
        return 1 if residue.get_id()[0] == " " else 0

    def accept_atom(self, atom) -> int:
        return 0 if atom.element in ("H", "D") else 1
