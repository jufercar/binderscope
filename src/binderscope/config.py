"""Configuration schema.

Everything target-specific lives in a YAML file, so the code itself carries no
knowledge of any particular target. A minimal run needs only ``designs`` and
``chains``; the ``reference`` section is optional and enables target
reconstruction (mutation reversion and reference-frame geometry).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

#: See binderscope.reversion for what each one does.
REVERSION_STRATEGIES = frozenset({"truncate", "donor"})


class ConfigError(ValueError):
    """Raised when a configuration file is structurally invalid."""


@dataclass
class Segment:
    """A contiguous stretch of the design-numbered target that maps onto a reference.

    ``input`` is the inclusive residue range as numbered in the design PDB;
    ``reference_start`` is the number the first residue of the range takes in
    the reference structure. The offset is derived, never written by hand.
    """

    input_start: int
    input_end: int
    reference_start: int
    label: str = ""

    @property
    def offset(self) -> int:
        return self.reference_start - self.input_start

    @property
    def length(self) -> int:
        return self.input_end - self.input_start + 1

    def __post_init__(self) -> None:
        if self.input_end < self.input_start:
            raise ConfigError(
                f"segment {self.label or '?'}: input_end < input_start "
                f"({self.input_end} < {self.input_start})"
            )


@dataclass
class Mutation:
    """A point substitution in the design-numbered target to revert to wild type."""

    position: int
    wild_type: str
    reference_position: int

    def __post_init__(self) -> None:
        self.wild_type = self.wild_type.upper()
        if self.wild_type not in AA3:
            raise ConfigError(
                f"mutation at {self.position}: '{self.wild_type}' is not a "
                "three-letter amino acid code"
            )


@dataclass
class ClashTarget:
    """A reference chain to measure binder clashes and volume overlap against."""

    chain: str
    label: str
    expect: str = "low"  # "low" | "high" — documentation only, drives no logic
    measure_volume: bool = False


@dataclass
class Reference:
    """Optional reference structure (e.g. a crystal) and its mapping."""

    structure: Path
    target_chain: str
    segments: list[Segment] = field(default_factory=list)
    linkers: list[tuple[int, int]] = field(default_factory=list)
    revert_mutations: list[Mutation] = field(default_factory=list)
    clash_targets: list[ClashTarget] = field(default_factory=list)
    native_target: Path | None = None

    @property
    def linker_positions(self) -> set[int]:
        out: set[int] = set()
        for lo, hi in self.linkers:
            out.update(range(lo, hi + 1))
        return out


@dataclass
class RosettaOptions:
    dalphaball: Path | None = None
    dssp: Path | None = None
    relax_repeats: int = 1
    relax_max_iter: int = 100
    relax_backbone: bool = False
    score_function: str = "beta_nov16"
    #: Interface specification for InterfaceAnalyzerMover. Left unset it is
    #: derived from the configured chains, which is what keeps the two from
    #: silently disagreeing.
    interface: str | None = None
    extra_flags: list[str] = field(default_factory=list)


@dataclass
class GeometryOptions:
    clash_hard: float = 2.0
    clash_soft: float = 2.5
    interface_cutoff: float = 4.0
    prerelax_clash_cutoff: float = 2.4
    atom_radius: float = 1.8
    probe_radius: float = 1.4
    voxel: float = 1.0


@dataclass
class RankMetric:
    col: str
    label: str
    higher_is_better: bool
    weight: float = 1.0
    description: str = ""


@dataclass
class Filter:
    col: str
    label: str
    op: str
    value: float
    enabled: bool = True

    def __post_init__(self) -> None:
        if self.op not in ("<", ">", "="):
            raise ConfigError(f"filter on '{self.col}': unsupported operator '{self.op}'")


@dataclass
class Config:
    name: str
    designs_dir: Path
    output_dir: Path
    target_chain: str = "A"
    binder_chain: str = "B"
    bindcraft_stats: Path | None = None
    reference: Reference | None = None
    rosetta: RosettaOptions = field(default_factory=RosettaOptions)
    geometry: GeometryOptions = field(default_factory=GeometryOptions)
    rank_metrics: list[RankMetric] = field(default_factory=list)
    filters: list[Filter] = field(default_factory=list)
    rank_exclude: list[Filter] = field(default_factory=list)
    pymol: dict[str, Any] = field(default_factory=dict)
    reversion_strategy: str = "truncate"
    description: str = ""

    # ── derived paths ────────────────────────────────────────────────────────
    @property
    def reverted_dir(self) -> Path:
        return self.output_dir / "reverted"

    @property
    def sessions_dir(self) -> Path:
        return self.output_dir / "sessions"

    @property
    def metrics_csv(self) -> Path:
        return self.output_dir / f"{self.name}_metrics.csv"

    @property
    def structural_csv(self) -> Path:
        return self.output_dir / f"{self.name}_metrics_structural.csv"

    @property
    def ranked_csv(self) -> Path:
        return self.output_dir / f"{self.name}_ranked.csv"

    @property
    def dashboard_html(self) -> Path:
        return self.output_dir / f"{self.name}_dashboard.html"

    @property
    def needs_reversion(self) -> bool:
        return bool(self.reference and self.reference.revert_mutations)

    @property
    def total_weight(self) -> float:
        return sum(m.weight for m in self.rank_metrics)


def _path(base: Path, value: Any) -> Path:
    """Resolve a config path relative to the config file's own directory."""
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def load_config(path: str | Path) -> Config:
    """Read and validate a YAML configuration file."""
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise ConfigError(f"config file not found: {path}")
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    if not isinstance(raw, dict):
        raise ConfigError(f"{path}: top level must be a mapping")
    return _build(raw, path.parent)


def _build(raw: dict[str, Any], base: Path) -> Config:
    for key in ("name", "designs", "output"):
        if key not in raw:
            raise ConfigError(f"missing required key: '{key}'")

    chains = raw.get("chains") or {}
    ref = _build_reference(raw.get("reference"), base)

    cfg = Config(
        name=str(raw["name"]),
        description=str(raw.get("description", "")),
        designs_dir=_path(base, raw["designs"]),
        output_dir=_path(base, raw["output"]),
        target_chain=str(chains.get("target", "A")),
        binder_chain=str(chains.get("binder", "B")),
        bindcraft_stats=(
            _path(base, raw["bindcraft_stats"]) if raw.get("bindcraft_stats") else None
        ),
        reference=ref,
        rosetta=_build_rosetta(raw.get("rosetta"), base),
        geometry=GeometryOptions(**(raw.get("geometry") or {})),
        rank_metrics=[
            RankMetric(
                col=m["col"],
                label=m.get("label", m["col"]),
                higher_is_better=bool(m["higher_is_better"]),
                weight=float(m.get("weight", 1.0)),
                description=m.get("description", ""),
            )
            for m in (raw.get("ranking") or {}).get("metrics", [])
        ],
        rank_exclude=[_build_filter(f) for f in (raw.get("ranking") or {}).get("exclude", [])],
        filters=[_build_filter(f) for f in raw.get("filters", [])],
        pymol=raw.get("pymol") or {},
        reversion_strategy=str(raw.get("reversion_strategy", "truncate")),
    )

    if cfg.needs_reversion and cfg.reference and not cfg.reference.native_target:
        raise ConfigError(
            "reference.revert_mutations requires reference.native_target "
            "(the as-designed target structure whose residues get reverted)"
        )

    if cfg.reversion_strategy not in REVERSION_STRATEGIES:
        raise ConfigError(
            f"reversion_strategy must be one of "
            f"{', '.join(sorted(REVERSION_STRATEGIES))}; got "
            f"{cfg.reversion_strategy!r}"
        )

    if "reversion_strategy" in (raw.get("pymol") or {}):
        raise ConfigError(
            "reversion_strategy belongs at the top level of the config, not "
            "under 'pymol' (it has nothing to do with session rendering)"
        )

    # The chains and the Rosetta interface specification describe the same two
    # chains. Letting them disagree would score an interface that is not there.
    derived = f"{cfg.target_chain}_{cfg.binder_chain}"
    if cfg.rosetta.interface is None:
        cfg.rosetta.interface = derived
    elif cfg.rosetta.interface != derived:
        raise ConfigError(
            f"rosetta.interface is '{cfg.rosetta.interface}' but chains are "
            f"target '{cfg.target_chain}' and binder '{cfg.binder_chain}', "
            f"which give '{derived}'. Remove rosetta.interface to derive it, "
            "or correct one of the two."
        )
    return cfg


def _build_filter(f: dict[str, Any]) -> Filter:
    return Filter(
        col=f["col"],
        label=f.get("label", f["col"]),
        op=str(f["op"]),
        value=float(f["value"]),
        enabled=bool(f.get("enabled", True)),
    )


def _build_rosetta(raw: dict[str, Any] | None, base: Path) -> RosettaOptions:
    raw = dict(raw or {})
    for key in ("dalphaball", "dssp"):
        if raw.get(key):
            raw[key] = _path(base, raw[key])
        else:
            raw.pop(key, None)
    return RosettaOptions(**raw)


def _build_reference(raw: dict[str, Any] | None, base: Path) -> Reference | None:
    if not raw:
        return None
    if "structure" not in raw or "target_chain" not in raw:
        raise ConfigError("reference requires 'structure' and 'target_chain'")

    segments = [
        Segment(
            input_start=int(s["input"][0]),
            input_end=int(s["input"][1]),
            reference_start=int(s["reference_start"]),
            label=s.get("label", f"seg{i + 1}"),
        )
        for i, s in enumerate(raw.get("segments", []))
    ]
    mutations = [
        Mutation(
            position=int(m["position"]),
            wild_type=str(m["wild_type"]),
            reference_position=int(m["reference_position"]),
        )
        for m in raw.get("revert_mutations", [])
    ]
    clash_targets = [
        ClashTarget(
            chain=str(c["chain"]),
            label=c.get("label", f"chain{c['chain']}"),
            expect=c.get("expect", "low"),
            measure_volume=bool(c.get("measure_volume", False)),
        )
        for c in raw.get("clash_targets", [])
    ]
    return Reference(
        structure=_path(base, raw["structure"]),
        target_chain=str(raw["target_chain"]),
        segments=segments,
        linkers=[(int(a), int(b)) for a, b in raw.get("linkers", [])],
        revert_mutations=mutations,
        clash_targets=clash_targets,
        native_target=(
            _path(base, raw["native_target"]) if raw.get("native_target") else None
        ),
    )
