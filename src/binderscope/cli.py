"""Command-line interface.

Stages are separate commands because they have very different costs: scoring is
hours of PyRosetta, while ranking and the dashboard are seconds and get re-run
many times as thresholds are tuned.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .config import ConfigError, load_config
from .dashboard import build_dashboard
from .pipeline import load_metrics, run_merge, run_reversion, run_scoring
from .ranking import apply_filters, rank_designs
from .sequences import extract_sequences


def _add_config_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("config", type=Path, help="path to the YAML run configuration")


def build_parser() -> argparse.ArgumentParser:
    """Assemble the argument parser, one subcommand per pipeline stage."""
    parser = argparse.ArgumentParser(
        prog="binderscope",
        description="Structure-aware triage of protein binder designs.",
    )
    parser.add_argument("--version", action="version", version=f"binderscope {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    revert = subparsers.add_parser(
        "revert", help="rebuild the wild-type target for every design"
    )
    _add_config_argument(revert)
    revert.add_argument(
        "--sessions", action="store_true", help="also write one PyMOL session per design"
    )

    score = subparsers.add_parser(
        "score", help="score interfaces with PyRosetta and measure placement"
    )
    _add_config_argument(score)

    merge = subparsers.add_parser(
        "merge", help="join BindCraft predictor statistics onto the metrics table"
    )
    _add_config_argument(merge)

    rank = subparsers.add_parser("rank", help="rank designs on the configured criteria")
    _add_config_argument(rank)
    rank.add_argument("--top", type=int, default=10, help="rows to print (default: 10)")
    rank.add_argument(
        "--metrics",
        type=Path,
        help="read this metrics CSV instead of the one implied by the config",
    )
    rank.add_argument(
        "--contributions",
        action="store_true",
        help="also print the per-metric contribution breakdown",
    )

    dashboard = subparsers.add_parser("dashboard", help="render the interactive dashboard")
    _add_config_argument(dashboard)
    dashboard.add_argument("--metrics", type=Path, help="metrics CSV to visualise")
    dashboard.add_argument("-o", "--output", type=Path, help="output HTML path")

    run = subparsers.add_parser(
        "run", help="revert (if configured), score, merge, rank and build the dashboard"
    )
    _add_config_argument(run)
    run.add_argument("--sessions", action="store_true", help="also write PyMOL sessions")

    sequences = subparsers.add_parser(
        "sequences", help="extract chain sequences from a directory of PDBs"
    )
    sequences.add_argument("pdb_dir", type=Path)
    sequences.add_argument("-c", "--chain", help="restrict to one chain, e.g. B")
    sequences.add_argument("-o", "--output", type=Path, help="FASTA output path")
    sequences.add_argument("--csv", type=Path, help="CSV output path")

    return parser


# ── commands ─────────────────────────────────────────────────────────────────
def cmd_revert(args) -> int:
    config = load_config(args.config)
    print(f"Reverting target for designs in {config.designs_dir}")
    report = run_reversion(config, sessions=args.sessions)
    print(f"Reversion: {report.summary()}")
    print(f"  structures: {config.reverted_dir}")
    return 1 if report.failures and not report.processed else 0


def cmd_score(args) -> int:
    config = load_config(args.config)
    print("Scoring interfaces (this is the slow stage; it resumes if interrupted)")
    report = run_scoring(config)
    print(f"Scoring: {report.summary()}")
    print(f"  metrics: {config.structural_csv}")
    return 1 if report.failures and not report.processed else 0


def cmd_merge(args) -> int:
    config = load_config(args.config)
    merged = run_merge(config)
    print(f"Merged table: {len(merged)} rows × {len(merged.columns)} columns")
    print(f"  written to: {config.metrics_csv}")
    return 0


def cmd_rank(args) -> int:
    config = load_config(args.config)
    df = pd.read_csv(args.metrics) if args.metrics else load_metrics(config)

    if config.filters:
        passing = apply_filters(df, config.filters)
        print(f"Designs passing all enabled filters: {passing.sum()} / {len(df)}")

    ranked, contributions = rank_designs(df, config.rank_metrics, config.rank_exclude)
    if ranked.empty:
        print("No designs left after the ranking exclusions.")
        return 1

    missing = ranked.attrs.get("missing_metrics")
    if missing:
        print(f"Note: configured metrics absent from the data: {', '.join(missing)}")

    ranked.to_csv(config.ranked_csv, index=False)

    columns = ["rank", "design", "score"] + [
        m.col for m in config.rank_metrics if m.col in ranked.columns
    ]
    top = ranked.head(args.top)
    print(f"\nTop {len(top)} of {len(ranked)} (max possible score {config.total_weight:g})")
    print(top[columns].round(3).to_string(index=False))

    if args.contributions:
        breakdown = contributions.head(args.top).round(3)
        breakdown.index = top["design"].to_numpy()
        print("\nPer-metric contributions")
        print(breakdown.to_string())

    print(f"\nFull ranking: {config.ranked_csv}")
    return 0


def cmd_dashboard(args) -> int:
    config = load_config(args.config)
    df = pd.read_csv(args.metrics) if args.metrics else load_metrics(config)
    output = build_dashboard(df, config, args.output)
    size_kb = output.stat().st_size / 1024
    print(f"Dashboard: {output} ({size_kb:.0f} KB, {len(df)} designs)")
    return 0


def cmd_run(args) -> int:
    config = load_config(args.config)

    if config.needs_reversion:
        print("[1/4] reversion")
        print(f"  {run_reversion(config, sessions=args.sessions).summary()}")
    else:
        print("[1/4] reversion — not configured, skipping")

    print("[2/4] scoring")
    print(f"  {run_scoring(config).summary()}")

    print("[3/4] merge")
    merged = run_merge(config)
    print(f"  {len(merged)} rows × {len(merged.columns)} columns")

    print("[4/4] ranking and dashboard")
    ranked, _ = rank_designs(merged, config.rank_metrics, config.rank_exclude)
    ranked.to_csv(config.ranked_csv, index=False)
    output = build_dashboard(merged, config)
    print(f"  ranking:   {config.ranked_csv}")
    print(f"  dashboard: {output}")
    return 0


def cmd_sequences(args) -> int:
    fasta = args.output or (args.pdb_dir / "sequences.fasta")
    records = extract_sequences(
        args.pdb_dir, chain=args.chain, fasta_out=fasta, csv_out=args.csv
    )
    print(f"{len(records)} sequences written to {fasta}")
    if args.csv:
        print(f"Summary table: {args.csv}")
    return 0


COMMANDS = {
    "revert": cmd_revert,
    "score": cmd_score,
    "merge": cmd_merge,
    "rank": cmd_rank,
    "dashboard": cmd_dashboard,
    "run": cmd_run,
    "sequences": cmd_sequences,
}


def main(argv: list[str] | None = None) -> int:
    """Entry point. Returns a process exit status rather than raising.

    Expected problems (a bad config, a missing input) are reported as a single
    line on stderr; an interrupt leaves completed work on disk so the run can
    be resumed.
    """
    args = build_parser().parse_args(argv)
    try:
        return COMMANDS[args.command](args)
    except (ConfigError, FileNotFoundError, NotADirectoryError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\ninterrupted; completed work is saved and the run will resume", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
