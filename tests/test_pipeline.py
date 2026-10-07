"""Stage orchestration: discovery, resumability and the merge step.

The scoring stage itself needs PyRosetta and is verified by running it on real
complexes. Everything around it — which designs get picked up, what happens on
restart, how a failure is recorded — is testable with synthetic structures, and
is what makes a multi-hour run survivable.
"""

from __future__ import annotations

import pandas as pd
import pytest
import yaml
from conftest import BINDER, TARGET_ENGINEERED, TARGET_REFERENCE, write_pdb

from binderscope.config import load_config
from binderscope.pipeline import (
    StageReport,
    find_designs,
    load_metrics,
    run_merge,
    run_reversion,
)


def build_run(tmp_path, *, n_designs=3, with_reversion=True, extra=None):
    """Write a small but complete run directory and return its loaded config."""
    designs = tmp_path / "designs"
    for i in range(n_designs):
        write_pdb(designs / f"design{i}_model1.pdb", {"A": TARGET_ENGINEERED, "B": BINDER})

    raw = {
        "name": "unit",
        "designs": "./designs",
        "output": "./analysis",
        "chains": {"target": "A", "binder": "B"},
    }
    if with_reversion:
        write_pdb(tmp_path / "native.pdb", {"D": TARGET_ENGINEERED})
        write_pdb(tmp_path / "reference.pdb", {"D": TARGET_REFERENCE})
        raw["reference"] = {
            "structure": "./reference.pdb",
            "target_chain": "D",
            "native_target": "./native.pdb",
            "segments": [{"input": [1, 10], "reference_start": 101}],
            "revert_mutations": [
                {"position": 5, "wild_type": "ALA", "reference_position": 105}
            ],
        }
    raw.update(extra or {})
    path = tmp_path / "run.yaml"
    path.write_text(yaml.safe_dump(raw))
    return load_config(path)


# ── discovery ────────────────────────────────────────────────────────────────
def test_designs_are_discovered_in_a_stable_order(tmp_path):
    config = build_run(tmp_path, n_designs=4, with_reversion=False)
    found = find_designs(config)
    assert len(found) == 4
    assert found == sorted(found)


def test_a_missing_designs_directory_is_reported(tmp_path):
    config = build_run(tmp_path, with_reversion=False)
    for pdb in config.designs_dir.glob("*.pdb"):
        pdb.unlink()
    config.designs_dir.rmdir()
    with pytest.raises(NotADirectoryError, match="designs directory not found"):
        find_designs(config)


def test_an_empty_designs_directory_is_reported(tmp_path):
    config = build_run(tmp_path, with_reversion=False)
    for pdb in config.designs_dir.glob("*.pdb"):
        pdb.unlink()
    with pytest.raises(FileNotFoundError, match="no PDB files"):
        find_designs(config)


# ── reversion stage and resumability ─────────────────────────────────────────
def test_reversion_processes_every_design(tmp_path):
    config = build_run(tmp_path, n_designs=3)
    report = run_reversion(config, verbose=False)
    assert report.processed == 3
    assert report.skipped == 0
    assert not report.failures
    assert len(list(config.reverted_dir.glob("*.pdb"))) == 3


def test_reversion_records_what_it_did(tmp_path):
    config = build_run(tmp_path, n_designs=2)
    run_reversion(config, verbose=False)
    log = pd.read_csv(config.reverted_dir / "reversion.csv")
    assert len(log) == 2
    assert set(log.columns) == {
        "design",
        "rmsd_design_native",
        "binder_len",
        "n_anchors",
        "n_reverted",
    }
    assert (log["n_reverted"] == 1).all()
    assert (log["n_anchors"] == 10).all()


def test_rerunning_skips_completed_work(tmp_path):
    # A hundred designs is hours of work; a restart must not redo any of it.
    config = build_run(tmp_path, n_designs=3)
    assert run_reversion(config, verbose=False).processed == 3
    second = run_reversion(config, verbose=False)
    assert second.processed == 0
    assert second.skipped == 3


def test_only_new_designs_are_processed_after_a_restart(tmp_path):
    config = build_run(tmp_path, n_designs=2)
    run_reversion(config, verbose=False)
    write_pdb(config.designs_dir / "design9_model1.pdb", {"A": TARGET_ENGINEERED, "B": BINDER})
    report = run_reversion(config, verbose=False)
    assert report.processed == 1
    assert report.skipped == 2


def test_a_deleted_output_is_regenerated(tmp_path):
    # The log alone is not trusted: the structure has to be on disk too.
    config = build_run(tmp_path, n_designs=2)
    run_reversion(config, verbose=False)
    next(iter(sorted(config.reverted_dir.glob("*.pdb")))).unlink()
    assert run_reversion(config, verbose=False).processed == 1


def test_one_bad_design_does_not_end_the_run(tmp_path):
    config = build_run(tmp_path, n_designs=2)
    # A design with no binder chain cannot be reverted.
    write_pdb(config.designs_dir / "broken_model1.pdb", {"A": TARGET_ENGINEERED})
    report = run_reversion(config, verbose=False)
    assert report.processed == 2
    assert len(report.failures) == 1
    assert report.failures[0][0] == "broken_model1"
    assert "not found" in report.failures[0][1]


def test_reversion_refuses_a_config_that_defines_none(tmp_path):
    config = build_run(tmp_path, with_reversion=False)
    with pytest.raises(ValueError, match="no reference.revert_mutations"):
        run_reversion(config, verbose=False)


# ── merge and loading ────────────────────────────────────────────────────────
def write_structural(config, designs) -> None:
    config.structural_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "design": designs,
            "design_base": [d.replace("_model1", "") for d in designs],
            "dG": [-20.0 - i for i in range(len(designs))],
        }
    ).to_csv(config.structural_csv, index=False)


def test_merge_without_predictor_stats_just_copies(tmp_path):
    config = build_run(tmp_path, with_reversion=False)
    write_structural(config, ["a_model1", "b_model1"])
    merged = run_merge(config)
    assert len(merged) == 2
    assert config.metrics_csv.is_file()


def test_merge_attaches_predictor_columns(tmp_path):
    stats = tmp_path / "stats.csv"
    pd.DataFrame(
        {"Design": ["a", "b"], "Average_i_pTM": [0.8, 0.6], "Sequence": ["AA", "CC"]}
    ).to_csv(stats, index=False)
    config = build_run(
        tmp_path, with_reversion=False, extra={"bindcraft_stats": "./stats.csv"}
    )
    write_structural(config, ["a_model1", "b_model1"])
    merged = run_merge(config)
    assert list(merged["Average_i_pTM"]) == [0.8, 0.6]
    assert "Sequence" in merged.columns


def test_merge_before_scoring_is_reported(tmp_path):
    config = build_run(tmp_path, with_reversion=False)
    with pytest.raises(FileNotFoundError, match="run scoring first"):
        run_merge(config)


def test_merged_table_is_preferred_over_the_structural_one(tmp_path):
    config = build_run(tmp_path, with_reversion=False)
    write_structural(config, ["a_model1"])
    assert len(load_metrics(config)) == 1
    run_merge(config)
    config.structural_csv.unlink()
    assert len(load_metrics(config)) == 1


def test_loading_metrics_before_any_run_is_reported(tmp_path):
    config = build_run(tmp_path, with_reversion=False)
    with pytest.raises(FileNotFoundError, match="no metrics found"):
        load_metrics(config)


# ── reporting ────────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    ("report", "expected"),
    [
        (StageReport(processed=3), "3 processed"),
        (StageReport(processed=1, skipped=2), "1 processed, 2 already done"),
        (
            StageReport(processed=1, failures=[("x", "boom")]),
            "1 processed, 1 failed",
        ),
    ],
)
def test_stage_summary_reads_clearly(report, expected):
    assert report.summary() == expected
