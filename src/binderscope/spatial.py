"""Geometry of a design placed back into a reference frame.

Interface energetics say how good a complex is on its own terms. They say
nothing about whether the binder sits where it was meant to sit. Superposing the
design onto a reference structure and measuring clashes against chosen reference
chains answers that second question: a binder that penetrates its target is
wrong however well it scores, and a binder that occupies the same space as a
native partner is evidence it engages the intended surface.

Which chains to measure against is configuration, not code — see
``reference.clash_targets``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.PDB import PDBParser, Superimposer

from .config import ClashTarget, GeometryOptions, Reference
from .geometry import count_contacting_atoms, heavy_atoms, volume_overlap
from .mapping import SegmentMap


class SpatialError(RuntimeError):
    """Raised when a design cannot be placed in the reference frame."""


class ReferenceFrame:
    """Loads a reference structure once and places designs into its frame."""

    def __init__(
        self,
        reference: Reference,
        segment_map: SegmentMap,
        geometry: GeometryOptions | None = None,
    ):
        self.reference = reference
        self.segment_map = segment_map
        self.geometry = geometry or GeometryOptions()

        model = PDBParser(QUIET=True).get_structure("reference", str(reference.structure))[0]
        if reference.target_chain not in model:
            raise SpatialError(
                f"reference chain '{reference.target_chain}' not found in "
                f"{reference.structure}"
            )
        self._model = model
        self._target_index = {
            r.get_id()[1]: r
            for r in model[reference.target_chain]
            if r.get_id()[0] == " "
        }

        self._clash_coords: dict[str, Any] = {}
        for target in reference.clash_targets:
            if target.chain not in model:
                raise SpatialError(
                    f"clash target chain '{target.chain}' not found in "
                    f"{reference.structure}"
                )
            self._clash_coords[target.label] = heavy_atoms(model[target.chain])

    def measure(
        self,
        pdb_file: str | Path,
        target_chain: str = "A",
        binder_chain: str = "B",
    ) -> dict[str, float]:
        """Superpose a design onto the reference and measure binder placement."""
        design = PDBParser(QUIET=True).get_structure("design", str(pdb_file))[0]
        for chain_id in (target_chain, binder_chain):
            if chain_id not in design:
                raise SpatialError(f"chain '{chain_id}' not found in {pdb_file}")

        design_index = {
            r.get_id()[1]: r for r in design[target_chain] if r.get_id()[0] == " "
        }

        fixed, moving = [], []
        for design_pos, reference_pos in self.segment_map.pairs:
            design_res = design_index.get(design_pos)
            reference_res = self._target_index.get(reference_pos)
            if design_res is None or reference_res is None:
                continue
            if "CA" in design_res and "CA" in reference_res:
                fixed.append(reference_res["CA"])
                moving.append(design_res["CA"])

        if len(fixed) < 3:
            raise SpatialError(
                f"only {len(fixed)} Cα anchors between design and reference "
                "(need at least 3)"
            )

        superimposer = Superimposer()
        superimposer.set_atoms(fixed, moving)
        superimposer.apply(list(design.get_atoms()))

        binder_xyz = heavy_atoms(design[binder_chain])
        results: dict[str, float] = {
            "rmsd_reference": float(superimposer.rms),
            "n_superposition_anchors": len(fixed),
        }

        for target in self.reference.clash_targets:
            results.update(self._measure_one(target, binder_xyz))
        return results

    def _measure_one(self, target: ClashTarget, binder_xyz) -> dict[str, float]:
        other = self._clash_coords[target.label]
        out: dict[str, float] = {
            f"clash_{target.label}_hard": count_contacting_atoms(
                binder_xyz, other, self.geometry.clash_hard
            ),
            f"clash_{target.label}_soft": count_contacting_atoms(
                binder_xyz, other, self.geometry.clash_soft
            ),
        }
        if target.measure_volume:
            overlap = volume_overlap(
                binder_xyz,
                other,
                atom_radius=self.geometry.atom_radius,
                probe_radius=self.geometry.probe_radius,
                voxel=self.geometry.voxel,
            )
            # The binder's own volume does not depend on what it is compared
            # against, so it is reported once under a target-free name.
            out["V_binder_A3"] = overlap["volume"]
            out[f"V_inter_{target.label}_A3"] = overlap["intersection"]
            out[f"V_inter_{target.label}_pct"] = overlap["intersection_pct"]
        return out
