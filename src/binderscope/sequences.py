"""Extract binder sequences from design structures."""

from __future__ import annotations

import csv
import warnings
from pathlib import Path

from Bio.PDB import PDBParser, PPBuilder


def extract_sequences(
    pdb_dir: str | Path,
    chain: str | None = None,
    fasta_out: str | Path | None = None,
    csv_out: str | Path | None = None,
) -> list[dict[str, object]]:
    """Read every PDB in a directory and collect its per-chain sequences.

    ``chain`` restricts the output to one chain (typically the binder). Returns
    the records, and writes FASTA and/or CSV when given a path.

    Raises ``ValueError`` when nothing could be extracted, rather than writing
    an empty FASTA that looks like a successful run.
    """
    pdb_dir = Path(pdb_dir)
    if not pdb_dir.is_dir():
        raise NotADirectoryError(f"not a directory: {pdb_dir}")

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"no PDB files in {pdb_dir}")

    parser = PDBParser(QUIET=True)
    builder = PPBuilder()
    records: list[dict[str, object]] = []

    for pdb_file in pdb_files:
        name = pdb_file.stem
        try:
            structure = parser.get_structure(name, str(pdb_file))
        except Exception as exc:
            warnings.warn(f"skipped {name}: cannot parse ({exc})", stacklevel=2)
            continue

        for chain_obj in structure[0]:
            chain_id = chain_obj.get_id()
            if chain and chain_id != chain:
                continue
            sequence = "".join(
                str(peptide.get_sequence()) for peptide in builder.build_peptides(chain_obj)
            )
            if sequence:
                records.append(
                    {
                        "design": name,
                        "chain": chain_id,
                        "length": len(sequence),
                        "sequence": sequence,
                    }
                )

    if not records:
        where = f"chain '{chain}' of " if chain else ""
        raise ValueError(
            f"no sequences could be extracted from {where}{len(pdb_files)} PDB "
            f"file(s) in {pdb_dir}; check the chain id and that the files hold "
            "connected polypeptide backbones"
        )

    if fasta_out:
        fasta_out = Path(fasta_out)
        fasta_out.parent.mkdir(parents=True, exist_ok=True)
        with open(fasta_out, "w") as fh:
            for record in records:
                fh.write(f">{record['design']}_chain{record['chain']}\n{record['sequence']}\n")

    if csv_out:
        csv_out = Path(csv_out)
        csv_out.parent.mkdir(parents=True, exist_ok=True)
        with open(csv_out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["design", "chain", "length", "sequence"])
            writer.writeheader()
            writer.writerows(records)

    return records
