"""Synthetic structures for tests, so no experimental data is needed."""

from __future__ import annotations

from pathlib import Path

import pytest

#: Heavy atoms per residue type, enough for the reversion tests.
ATOMS = {
    "GLY": ["N", "CA", "C", "O"],
    "ALA": ["N", "CA", "C", "O", "CB"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
}


def residue_lines(
    serial: int,
    resname: str,
    chain: str,
    number: int,
    origin: tuple[float, float, float],
) -> tuple[list[str], int]:
    """PDB ATOM records for one residue, atoms fanned out from ``origin``."""
    lines = []
    x0, y0, z0 = origin
    for offset, name in enumerate(ATOMS[resname]):
        x, y, z = x0 + offset * 0.8, y0, z0
        element = name[0]
        # Columns matter: 13-16 atom name, 17 altLoc, 18-20 resName,
        # 22 chainID, 23-26 resSeq.
        lines.append(
            f"ATOM  {serial:5d} {name:<4s} {resname:>3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
        )
        serial += 1
    return lines, serial


def write_pdb(path: Path, chains: dict[str, list[tuple[int, str]]]) -> Path:
    """Write a PDB from ``{chain_id: [(residue_number, resname), ...]}``.

    Residues are spaced 3.8 Å apart along x and chains offset along y, so
    superposition has well-defined geometry to work with.
    """
    serial = 1
    lines: list[str] = []
    for chain_index, (chain_id, residues) in enumerate(chains.items()):
        for residue_index, (number, resname) in enumerate(residues):
            origin = (residue_index * 3.8, chain_index * 12.0, 0.0)
            new_lines, serial = residue_lines(serial, resname, chain_id, number, origin)
            lines.extend(new_lines)
        lines.append("TER\n")
    lines.append("END\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))
    return path


#: Ten target residues, position 5 engineered to lysine, plus a short binder.
TARGET_ENGINEERED = [
    (n, "LYS" if n == 5 else "GLY") for n in range(1, 11)
]
BINDER = [(n, "ALA") for n in range(1, 6)]
#: The same stretch in reference numbering, with the wild-type residue restored.
TARGET_REFERENCE = [
    (100 + n, "ALA" if n == 5 else "GLY") for n in range(1, 11)
]


@pytest.fixture
def design_pdb(tmp_path) -> Path:
    return write_pdb(tmp_path / "design.pdb", {"A": TARGET_ENGINEERED, "B": BINDER})


@pytest.fixture
def native_target_pdb(tmp_path) -> Path:
    return write_pdb(tmp_path / "native.pdb", {"D": TARGET_ENGINEERED})


@pytest.fixture
def reference_pdb(tmp_path) -> Path:
    return write_pdb(
        tmp_path / "reference.pdb",
        {"D": TARGET_REFERENCE, "B": [(200 + n, "VAL") for n in range(1, 6)]},
    )
