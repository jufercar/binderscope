"""Interface composition and secondary structure, from coordinates alone."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio.PDB import DSSP, PDBParser
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

from .geometry import HYDROGENS

HYDROPHOBIC = set("ACFGILMPVWY")


def interface_residues(
    pdb_file: str | Path,
    target_chain: str = "A",
    binder_chain: str = "B",
    cutoff: float = 4.0,
) -> dict[int, str]:
    """Binder residues with any heavy atom within ``cutoff`` of the target.

    Returns ``{residue_number: one_letter_code}``.
    """
    structure = PDBParser(QUIET=True).get_structure("tmp", str(pdb_file))
    model = structure[0]
    if target_chain not in model or binder_chain not in model:
        return {}

    target_coords = [
        atom.get_vector().get_array()
        for residue in model[target_chain]
        if residue.get_id()[0] == " "
        for atom in residue
        if atom.element not in HYDROGENS
    ]
    if not target_coords:
        return {}
    target_tree = cKDTree(np.asarray(target_coords, dtype=np.float32))

    found: dict[int, str] = {}
    for residue in model[binder_chain]:
        if residue.get_id()[0] != " ":
            continue
        coords = [
            atom.get_vector().get_array()
            for atom in residue
            if atom.element not in HYDROGENS
        ]
        if not coords:
            continue
        hits = target_tree.query_ball_point(
            np.asarray(coords, dtype=np.float32), r=cutoff, return_length=True
        )
        if np.any(hits > 0):
            found[residue.get_id()[1]] = seq1(residue.get_resname())
    return found


def _classify(code: str) -> str:
    if code in ("H", "G", "I"):
        return "H"
    if code == "E":
        return "E"
    return "L"


def _percentages(codes: list[str]) -> tuple[float, float, float]:
    n = len(codes) or 1
    return (
        codes.count("H") / n * 100,
        codes.count("E") / n * 100,
        codes.count("L") / n * 100,
    )


def secondary_structure(
    pdb_file: str | Path,
    binder_chain: str = "B",
    interface_positions: list[int] | None = None,
    dssp_binary: str | Path | None = None,
) -> dict[str, float | None]:
    """Helix/sheet/loop percentages for the whole binder and for its interface.

    Returns all-``None`` values when DSSP is unavailable or fails, so a missing
    binary degrades the metric table instead of aborting the run.
    """
    empty: dict[str, float | None] = {
        "binder_helix_pct": None,
        "binder_sheet_pct": None,
        "binder_loop_pct": None,
        "iface_helix_pct": None,
        "iface_sheet_pct": None,
        "iface_loop_pct": None,
    }

    structure = PDBParser(QUIET=True).get_structure("tmp", str(pdb_file))
    model = structure[0]
    try:
        dssp = (
            DSSP(model, str(pdb_file), dssp=str(dssp_binary))
            if dssp_binary
            else DSSP(model, str(pdb_file))
        )
    except Exception:
        return empty

    ss_map = {(key[0], key[1][1]): value[2] for key, value in dssp.property_dict.items()}
    if binder_chain not in model:
        return empty

    binder_codes = [
        _classify(ss_map.get((binder_chain, residue.get_id()[1]), "-"))
        for residue in model[binder_chain]
        if residue.get_id()[0] == " "
    ]
    helix, sheet, loop = _percentages(binder_codes)

    if interface_positions:
        iface_codes = [
            _classify(ss_map.get((binder_chain, pos), "-")) for pos in interface_positions
        ]
        i_helix, i_sheet, i_loop = _percentages(iface_codes)
    else:
        i_helix = i_sheet = i_loop = None

    return {
        "binder_helix_pct": helix,
        "binder_sheet_pct": sheet,
        "binder_loop_pct": loop,
        "iface_helix_pct": i_helix,
        "iface_sheet_pct": i_sheet,
        "iface_loop_pct": i_loop,
    }


def hydrophobic_fraction(residues: dict[int, str] | list[str]) -> float:
    """Fraction of the given residues whose one-letter code is hydrophobic."""
    codes = list(residues.values()) if isinstance(residues, dict) else list(residues)
    if not codes:
        return 0.0
    return sum(1 for code in codes if code in HYDROPHOBIC) / len(codes)
