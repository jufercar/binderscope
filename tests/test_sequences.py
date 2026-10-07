"""Sequence extraction.

Needs a backbone Biopython will actually chain together, so these build their
own structures with plausible peptide geometry rather than reusing the
coarse fixtures used for superposition tests.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from binderscope.sequences import extract_sequences

#: Offsets along x within one residue, and the spacing to the next, chosen so
#: the C-N distance between consecutive residues is ~1.33 Å and PPBuilder links
#: them into one polypeptide.
BACKBONE_OFFSETS = {"N": 0.00, "CA": 1.46, "C": 2.40, "O": 2.60}
RESIDUE_SPACING = 3.73


def write_peptide(path: Path, chains: dict[str, str]) -> Path:
    """Write a PDB of connected poly-residue chains from one-letter sequences."""
    three = {
        "A": "ALA", "G": "GLY", "V": "VAL", "L": "LEU",
        "S": "SER", "K": "LYS", "E": "GLU", "T": "THR",
    }
    serial = 1
    lines: list[str] = []
    for chain_index, (chain_id, sequence) in enumerate(chains.items()):
        for position, code in enumerate(sequence, start=1):
            resname = three[code]
            base = (position - 1) * RESIDUE_SPACING
            for name, offset in BACKBONE_OFFSETS.items():
                x = base + offset
                y = chain_index * 25.0
                lines.append(
                    f"ATOM  {serial:5d} {name:<4s} {resname:>3s} "
                    f"{chain_id}{position:4d}    "
                    f"{x:8.3f}{y:8.3f}{0.0:8.3f}  1.00  0.00          "
                    f"{name[0]:>2s}\n"
                )
                serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))
    return path


@pytest.fixture
def designs(tmp_path) -> Path:
    directory = tmp_path / "designs"
    write_peptide(directory / "d1.pdb", {"A": "GAVLSKET", "B": "AAGVL"})
    write_peptide(directory / "d2.pdb", {"A": "GAVLSKET", "B": "LLKSE"})
    return directory


def test_every_chain_is_extracted_by_default(designs):
    records = extract_sequences(designs)
    assert len(records) == 4
    assert {r["chain"] for r in records} == {"A", "B"}


def test_a_single_chain_can_be_selected(designs):
    records = extract_sequences(designs, chain="B")
    assert [r["sequence"] for r in records] == ["AAGVL", "LLKSE"]
    assert all(r["length"] == 5 for r in records)


def test_records_name_their_design(designs):
    records = extract_sequences(designs, chain="B")
    assert [r["design"] for r in records] == ["d1", "d2"]


def test_fasta_output_is_written(designs, tmp_path):
    fasta = tmp_path / "out.fasta"
    extract_sequences(designs, chain="B", fasta_out=fasta)
    lines = fasta.read_text().splitlines()
    assert lines[0] == ">d1_chainB"
    assert lines[1] == "AAGVL"
    assert len(lines) == 4


def test_csv_output_is_written(designs, tmp_path):
    out = tmp_path / "out.csv"
    extract_sequences(designs, chain="B", csv_out=out)
    with open(out) as fh:
        rows = list(csv.DictReader(fh))
    assert [r["sequence"] for r in rows] == ["AAGVL", "LLKSE"]
    assert rows[0]["length"] == "5"


def test_an_absent_chain_is_an_error_not_an_empty_file(designs, tmp_path):
    # It used to write a 0-byte FASTA and report success, which is the kind of
    # thing you only notice after ordering nothing.
    fasta = tmp_path / "empty.fasta"
    with pytest.raises(ValueError, match="no sequences could be extracted"):
        extract_sequences(designs, chain="Z", fasta_out=fasta)
    assert not fasta.exists()


def test_a_missing_directory_is_reported(tmp_path):
    with pytest.raises(NotADirectoryError, match="not a directory"):
        extract_sequences(tmp_path / "nope")


def test_a_directory_without_pdbs_is_reported(tmp_path):
    (tmp_path / "empty").mkdir()
    with pytest.raises(FileNotFoundError, match="no PDB files"):
        extract_sequences(tmp_path / "empty")


def test_an_unparseable_file_warns_and_is_skipped(designs):
    (designs / "broken.pdb").write_text("ATOM      1  CA  ALA A   1    not-coords\nEND\n")
    with pytest.warns(UserWarning, match="cannot parse"):
        records = extract_sequences(designs, chain="B")
    assert len(records) == 2
