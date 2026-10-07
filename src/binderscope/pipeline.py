"""Stage orchestration.

Scoring a hundred designs takes hours, mostly in FastRelax, and a run that loses
everything to one malformed PDB is unusable in practice. Every stage is
therefore resumable: results are appended to a CSV as they are produced, and a
restart skips whatever is already there. A design that fails is recorded and the
run continues.
"""

from __future__ import annotations

import gc
import time
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .config import Config
from .mapping import SegmentMap
from .merge import merge_predictor_stats, split_design_name
from .reversion import ReversionError, revert_design
from .spatial import ReferenceFrame


@dataclass
class StageReport:
    processed: int = 0
    skipped: int = 0
    failures: list[tuple[str, str]] = field(default_factory=list)

    def summary(self) -> str:
        parts = [f"{self.processed} processed"]
        if self.skipped:
            parts.append(f"{self.skipped} already done")
        if self.failures:
            parts.append(f"{len(self.failures)} failed")
        return ", ".join(parts)


def find_designs(config: Config) -> list[Path]:
    """All design PDBs to analyse, sorted for reproducible ordering."""
    if not config.designs_dir.is_dir():
        raise NotADirectoryError(f"designs directory not found: {config.designs_dir}")
    designs = sorted(config.designs_dir.glob("*.pdb"))
    if not designs:
        raise FileNotFoundError(f"no PDB files in {config.designs_dir}")
    return designs


def _done_designs(csv_path: Path) -> set[str]:
    if not csv_path.is_file():
        return set()
    try:
        existing = pd.read_csv(csv_path, usecols=["design"])
    except (ValueError, pd.errors.EmptyDataError):
        return set()
    return set(existing["design"].astype(str))


def _append_row(csv_path: Path, row: dict) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame([row])
    frame.to_csv(
        csv_path,
        mode="a" if csv_path.is_file() else "w",
        header=not csv_path.is_file(),
        index=False,
    )


# ── stage 1: reversion ───────────────────────────────────────────────────────
def run_reversion(
    config: Config,
    sessions: bool = False,
    verbose: bool = True,
) -> StageReport:
    """Rebuild the wild-type target for every design."""
    if not config.needs_reversion:
        raise ValueError(
            "this config defines no reference.revert_mutations; "
            "the reversion stage has nothing to do"
        )
    assert config.reference is not None

    segment_map = SegmentMap.from_reference(config.reference)
    report = StageReport()
    metadata_csv = config.reverted_dir / "reversion.csv"
    done = _done_designs(metadata_csv)

    pymol_binary = None
    if sessions:
        from .pymol_session import find_pymol, write_session

        pymol_binary = find_pymol(config)
        if not pymol_binary and verbose:
            print("  PyMOL not found; skipping session generation")

    for design_pdb in find_designs(config):
        name = design_pdb.stem
        output_pdb = config.reverted_dir / f"{name}.pdb"
        if name in done and output_pdb.is_file():
            report.skipped += 1
            continue

        started = time.time()
        try:
            result = revert_design(
                design_pdb=design_pdb,
                output_pdb=output_pdb,
                native_target_pdb=config.reference.native_target,
                reference_pdb=config.reference.structure,
                segment_map=segment_map,
                mutations=config.reference.revert_mutations,
                native_chain=config.reference.target_chain,
                reference_chain=config.reference.target_chain,
                design_target_chain=config.target_chain,
                design_binder_chain=config.binder_chain,
                strategy=config.pymol.get("reversion_strategy", "truncate"),
            )
        except (ReversionError, KeyError, ValueError) as exc:
            report.failures.append((name, str(exc)))
            if verbose:
                print(f"  FAIL {name}: {exc}")
            continue

        _append_row(
            metadata_csv,
            {
                "design": name,
                "rmsd_design_native": round(result.rmsd_ca, 4),
                "binder_len": result.binder_length,
                "n_anchors": result.anchors,
                "n_reverted": len(result.reverted),
            },
        )
        report.processed += 1
        if verbose:
            print(
                f"  ok   {name} | RMSD {result.rmsd_ca:.3f} Å | "
                f"binder {result.binder_length} res | {time.time() - started:.1f}s"
            )

        if sessions and pymol_binary:
            write_session(
                output_pdb,
                config.sessions_dir / f"{name}.pse",
                config,
                segment_map=segment_map,
                pymol_binary=pymol_binary,
            )

    return report


# ── stage 2: scoring ─────────────────────────────────────────────────────────
def _scoring_inputs(config: Config) -> Iterator[Path]:
    if config.needs_reversion:
        if not config.reverted_dir.is_dir():
            raise FileNotFoundError(
                f"no reverted structures in {config.reverted_dir}; "
                "run the reversion stage first"
            )
        yield from sorted(config.reverted_dir.glob("*.pdb"))
    else:
        yield from find_designs(config)


def run_scoring(config: Config, verbose: bool = True) -> StageReport:
    """Score every complex with PyRosetta, plus placement geometry if configured."""
    from .rosetta import build_scorer

    report = StageReport()
    done = _done_designs(config.structural_csv)

    frame = None
    if config.reference and config.reference.clash_targets:
        frame = ReferenceFrame(
            config.reference,
            SegmentMap.from_reference(config.reference),
            config.geometry,
        )

    scorer = build_scorer(config)
    if verbose:
        print(
            "  BuriedUnsatHbonds: "
            + ("DAlphaBall" if scorer.use_dalphaball else "InterfaceAnalyzer fallback")
        )

    reversion_csv = config.reverted_dir / "reversion.csv"
    reversion_data: dict[str, dict] = {}
    if reversion_csv.is_file():
        reversion_data = (
            pd.read_csv(reversion_csv).set_index("design").to_dict(orient="index")
        )

    for pdb_file in _scoring_inputs(config):
        name = pdb_file.stem
        if name in done:
            report.skipped += 1
            continue

        base, model = split_design_name(name)
        started = time.time()
        row: dict = {"design": name, "design_base": base, "model": model}
        row.update(reversion_data.get(name, {}))

        try:
            row.update(scorer.score(pdb_file, target_chain=config.target_chain))
            if frame is not None:
                row.update(
                    frame.measure(
                        pdb_file,
                        target_chain=config.target_chain,
                        binder_chain=config.binder_chain,
                    )
                )
        except Exception as exc:
            report.failures.append((name, str(exc)))
            if verbose:
                print(f"  FAIL {name}: {exc}")
            continue

        _append_row(config.structural_csv, row)
        report.processed += 1
        if verbose:
            print(
                f"  ok   {name} | dG {row.get('dG', float('nan')):.1f} | "
                f"{time.time() - started:.0f}s"
            )
        gc.collect()

    return report


# ── stage 3: merge ───────────────────────────────────────────────────────────
def run_merge(config: Config) -> pd.DataFrame:
    """Join predictor statistics onto the structural metrics and write the table."""
    if not config.structural_csv.is_file():
        raise FileNotFoundError(
            f"no structural metrics at {config.structural_csv}; run scoring first"
        )
    metrics = pd.read_csv(config.structural_csv)
    merged = merge_predictor_stats(metrics, config.bindcraft_stats)
    merged.to_csv(config.metrics_csv, index=False)
    return merged


def load_metrics(config: Config) -> pd.DataFrame:
    """Best available metrics table: merged if present, else structural only."""
    for path in (config.metrics_csv, config.structural_csv):
        if path.is_file():
            return pd.read_csv(path)
    raise FileNotFoundError(
        f"no metrics found; expected {config.metrics_csv} or {config.structural_csv}"
    )
