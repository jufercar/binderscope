"""PyRosetta interface scoring.

PyRosetta is imported lazily so the rest of the package (ranking, dashboard,
geometry) stays importable and testable without it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.PDB import PDBParser

from .config import Config, GeometryOptions, RosettaOptions
from .geometry import count_interchain_clashes
from .interface import (
    HYDROPHOBIC,
    hydrophobic_fraction,
    interface_residues,
    secondary_structure,
)

_BUNS_XML = """
<ROSETTASCRIPTS>
 <SCOREFXNS><ScoreFunction name="scorefxn" weights="{weights}"/></SCOREFXNS>
 <FILTERS>
  <BuriedUnsatHbonds name="buns"
    report_all_heavy_atom_unsats="true"
    scorefxn="scorefxn"
    ignore_surface_res="false"
    use_ddG_style="true"
    dalphaball_sasa="1"
    probe_radius="1.1"
    burial_cutoff_apo="0.2"
    confidence="0" />
 </FILTERS>
 <PROTOCOLS/>
</ROSETTASCRIPTS>
"""


class RosettaUnavailable(RuntimeError):
    """Raised when PyRosetta is required but not importable."""


def _scalar(value: Any) -> float:
    """Extract a float from a PyRosetta scalar or vector1_double."""
    try:
        return float(value)
    except TypeError:
        return float(list(value)[0])


class InterfaceScorer:
    """Holds an initialised PyRosetta session and scores one complex at a time.

    Initialising PyRosetta is a process-global, one-shot operation, so a single
    scorer instance is reused across all designs in a run.
    """

    def __init__(
        self,
        options: RosettaOptions,
        geometry: GeometryOptions | None = None,
        binder_chain: str = "B",
    ):
        try:
            import pyrosetta
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise RosettaUnavailable(
                "PyRosetta is required for interface scoring. Install it, or run "
                "only the 'rank' and 'dashboard' stages on an existing metrics CSV."
            ) from exc

        self._pyrosetta = pyrosetta
        self.options = options
        self.geometry = geometry or GeometryOptions()
        self.binder_chain = binder_chain

        beta = "true" if options.score_function == "beta_nov16" else "false"
        flags = [
            f"-corrections:beta_nov16 {beta}",
            "-ignore_unrecognized_res",
            "-ignore_zero_occupancy",
            f"-relax:default_repeats {options.relax_repeats}",
            "-mute all",
        ]
        if options.dalphaball:
            flags.append(f"-holes:dalphaball {options.dalphaball}")
        flags.extend(options.extra_flags)
        pyrosetta.init(" ".join(flags))

        self.score_function = pyrosetta.get_fa_scorefxn()
        self.buns_unavailable_reason: str | None = None
        self._relax = self._build_relax()
        self._buns_filter, self.use_dalphaball = self._build_buns()

    # ── setup helpers ────────────────────────────────────────────────────────
    def _build_relax(self):
        from pyrosetta.rosetta.core.kinematics import MoveMap
        from pyrosetta.rosetta.protocols.relax import FastRelax

        relax = FastRelax(self.score_function, standard_repeats=self.options.relax_repeats)
        movemap = MoveMap()
        movemap.set_bb(self.options.relax_backbone)
        movemap.set_chi(True)
        movemap.set_jump(True)
        relax.set_movemap(movemap)
        relax.max_iter(self.options.relax_max_iter)
        return relax

    def _build_buns(self):
        """Build the BuriedUnsatHbonds filter, or explain why it is unavailable.

        The fallback used to be silent, which hid a wrong call signature behind
        plausible-looking numbers from a different estimator. Any failure is
        now recorded in ``buns_unavailable_reason`` and reported by the caller.
        """
        from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

        if not self.options.dalphaball:
            self.buns_unavailable_reason = "no dalphaball binary configured"
            return None, False

        if not Path(self.options.dalphaball).is_file():
            self.buns_unavailable_reason = (
                f"dalphaball binary not found at {self.options.dalphaball}"
            )
            return None, False

        try:
            xml = _BUNS_XML.format(weights=self.options.score_function)
            # Note: static_get_filter() takes the XML text alone — passing a
            # filter name as a second argument raises TypeError. Going through
            # XmlObjects keeps the configured score function meaningful.
            return XmlObjects.create_from_string(xml).get_filter("buns"), True
        except Exception as exc:
            self.buns_unavailable_reason = f"{type(exc).__name__}: {exc}"
            return None, False

    # ── scoring ──────────────────────────────────────────────────────────────
    def score(self, pdb_file: str | Path, target_chain: str = "A") -> dict[str, Any]:
        """Relax and score one complex, returning a flat metric dictionary."""
        from pyrosetta.rosetta.core.select.residue_selector import (
            ChainSelector,
            LayerSelector,
        )
        from pyrosetta.rosetta.core.simple_metrics import metrics as pr_metrics
        from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

        pdb_file = str(pdb_file)
        model = PDBParser(QUIET=True).get_structure("tmp", pdb_file)[0]
        if self.binder_chain not in model:
            raise ValueError(
                f"binder chain '{self.binder_chain}' not found in {pdb_file}"
            )
        binder_len = sum(
            1 for r in model[self.binder_chain] if r.get_id()[0] == " "
        )
        clashes_pre_relax = count_interchain_clashes(
            model, threshold=self.geometry.prerelax_clash_cutoff
        )

        pose = self._pyrosetta.pose_from_pdb(pdb_file)
        self._relax.apply(pose)

        analyzer = InterfaceAnalyzerMover(self.options.interface)
        analyzer.set_compute_packstat(True)
        analyzer.set_compute_interface_energy(True)
        analyzer.set_calc_dSASA(True)
        analyzer.set_calc_hbond_sasaE(True)
        analyzer.set_compute_interface_sc(True)
        analyzer.set_pack_separated(True)
        analyzer.apply(pose)
        data = analyzer.get_all_data()

        dG = _scalar(analyzer.get_interface_dG())
        dSASA = _scalar(data.dSASA)
        n_hbonds = _scalar(data.interface_hbonds)

        if self.use_dalphaball and self._buns_filter is not None:
            try:
                unsat = self._buns_filter.report_sm(pose)
            except Exception:
                unsat = analyzer.get_interface_delta_hbond_unsat()
        else:
            unsat = analyzer.get_interface_delta_hbond_unsat()

        iface = interface_residues(
            pdb_file,
            target_chain=target_chain,
            binder_chain=self.binder_chain,
            cutoff=self.geometry.interface_cutoff,
        )
        n_iface = len(iface) or 1

        selector = ChainSelector(self.binder_chain)
        energy_metric = pr_metrics.TotalEnergyMetric()
        energy_metric.set_residue_selector(selector)
        binder_score = energy_metric.calculate(pose)

        sasa_metric = pr_metrics.SasaMetric()
        sasa_metric.set_residue_selector(selector)
        binder_sasa = sasa_metric.calculate(pose)

        binder_pose = pose.split_by_chain(
            self._chain_number(pose, self.binder_chain)
        )
        layers = LayerSelector()
        layers.set_layers(False, False, True)  # surface only
        surface_mask = layers.apply(binder_pose)
        surface = [
            binder_pose.residue(i).name1()
            for i in range(1, binder_pose.total_residue() + 1)
            if surface_mask[i]
        ]

        ss = secondary_structure(
            pdb_file,
            binder_chain=self.binder_chain,
            interface_positions=list(iface),
            dssp_binary=self.options.dssp,
        )

        metrics: dict[str, Any] = {
            # Always reported, so a run without a reversion stage still gets it.
            "binder_len": binder_len,
            "clashes_pre_relax": clashes_pre_relax,
            "dG": dG,
            "dSASA": dSASA,
            "dG_dSASA_ratio": _scalar(data.dG_dSASA_ratio) * 100,
            "sc": _scalar(data.sc_value),
            "packstat": _scalar(data.packstat),
            "n_hbonds": n_hbonds,
            "hbond_pct": n_hbonds / n_iface * 100,
            "unsat_hbonds": unsat,
            "unsat_hbonds_pct": unsat / n_iface * 100,
            "n_interface_res": len(iface),
            "interface_hydro": hydrophobic_fraction(iface) * 100,
            "interface_AA": dict(iface),
            "binder_score": binder_score,
            "binder_sasa": binder_sasa,
            "interface_sasa_pct": dSASA / binder_sasa * 100 if binder_sasa else 0.0,
            "surface_hydro": (
                sum(1 for aa in surface if aa in HYDROPHOBIC) / len(surface)
                if surface
                else 0.0
            ),
        }
        metrics.update(ss)
        return metrics

    @staticmethod
    def _chain_number(pose, chain_id: str) -> int:
        """Rosetta chain index for a PDB chain letter (1-based)."""
        info = pose.pdb_info()
        for i in range(1, pose.num_chains() + 1):
            begin = pose.chain_begin(i)
            if info.chain(begin) == chain_id:
                return i
        raise ValueError(f"chain '{chain_id}' not present in pose")


def build_scorer(config: Config) -> InterfaceScorer:
    return InterfaceScorer(
        config.rosetta,
        geometry=config.geometry,
        binder_chain=config.binder_chain,
    )
