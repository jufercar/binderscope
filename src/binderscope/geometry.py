"""Coordinate-level geometry: heavy-atom extraction, clash counting, volume overlap.

All functions here are target-agnostic and take plain coordinate arrays, so they
are cheap to unit-test without any structure files.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

HYDROGENS = ("H", "D")


def heavy_atoms(chain, res_min: int | None = None, res_max: int | None = None) -> np.ndarray:
    """Heavy-atom coordinates of a Bio.PDB chain, optionally restricted to a range.

    Hydrogens and non-standard records (waters, ligands) are skipped. Returns an
    ``(n, 3)`` float32 array, empty if nothing matched.
    """
    coords: list[np.ndarray] = []
    for residue in chain:
        if residue.get_id()[0] != " ":
            continue
        number = residue.get_id()[1]
        if res_min is not None and number < res_min:
            continue
        if res_max is not None and number > res_max:
            continue
        for atom in residue:
            if atom.element not in HYDROGENS:
                coords.append(atom.get_vector().get_array())
    if not coords:
        return np.zeros((0, 3), dtype=np.float32)
    return np.asarray(coords, dtype=np.float32)


def count_contacting_atoms(xyz_a: np.ndarray, xyz_b: np.ndarray, threshold: float) -> int:
    """Number of atoms in ``xyz_a`` within ``threshold`` of any atom in ``xyz_b``.

    This counts *atoms of a*, not pairs, so the value is bounded by ``len(xyz_a)``
    and does not blow up when two dense regions interpenetrate.
    """
    if len(xyz_a) == 0 or len(xyz_b) == 0:
        return 0
    tree = cKDTree(xyz_b)
    counts = tree.query_ball_point(xyz_a, r=threshold, return_length=True)
    return int(np.count_nonzero(counts))


def count_interchain_clashes(structure_model, threshold: float = 2.4) -> int:
    """Heavy-atom pairs closer than ``threshold`` that belong to different chains."""
    coords: list[np.ndarray] = []
    chains: list[str] = []
    for chain in structure_model:
        for residue in chain:
            if residue.get_id()[0] != " ":
                continue
            for atom in residue:
                if atom.element not in HYDROGENS:
                    coords.append(atom.get_vector().get_array())
                    chains.append(chain.id)
    if not coords:
        return 0
    tree = cKDTree(np.asarray(coords, dtype=np.float32))
    return sum(1 for i, j in tree.query_pairs(threshold) if chains[i] != chains[j])


def occupancy_grid(
    coords: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    cut: float,
    voxel: float = 1.0,
) -> np.ndarray:
    """Boolean grid over ``[lower, upper]`` marking voxels within ``cut`` of an atom."""
    dims = [max(int(round((hi - lo) / voxel)) + 1, 1) for lo, hi in zip(lower, upper, strict=False)]
    if len(coords) == 0:
        return np.zeros(dims, dtype=bool)
    axes = [
        np.arange(n, dtype=np.float32) * voxel + lo for n, lo in zip(dims, lower, strict=False)
    ]
    grid = np.meshgrid(*axes, indexing="ij")
    points = np.column_stack([g.ravel() for g in grid])
    tree = cKDTree(coords)
    hits = tree.query_ball_point(points, r=cut, return_length=True)
    return (hits > 0).reshape(dims)


def volume_overlap(
    xyz_query: np.ndarray,
    xyz_other: np.ndarray,
    atom_radius: float = 1.8,
    probe_radius: float = 1.4,
    voxel: float = 1.0,
) -> dict[str, float]:
    """Voxelised volume of ``xyz_query`` and its intersection with ``xyz_other``.

    The grid is bounded by the query's own extent plus padding, so the reported
    percentage is always "fraction of the query volume that is also occupied by
    the other selection".
    """
    if len(xyz_query) == 0:
        return {"volume": 0.0, "intersection": 0.0, "intersection_pct": 0.0}

    cut = atom_radius + probe_radius
    pad = cut + 3.0
    lower = xyz_query.min(axis=0) - pad
    upper = xyz_query.max(axis=0) + pad

    grid_query = occupancy_grid(xyz_query, lower, upper, cut, voxel)
    grid_other = occupancy_grid(xyz_other, lower, upper, cut, voxel)

    unit = voxel**3
    volume = float(grid_query.sum()) * unit
    intersection = float((grid_query & grid_other).sum()) * unit
    pct = 100.0 * intersection / volume if volume > 0 else 0.0
    return {"volume": volume, "intersection": intersection, "intersection_pct": pct}
