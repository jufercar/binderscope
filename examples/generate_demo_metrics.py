#!/usr/bin/env python3
"""Generate the synthetic metrics table used by the demo.

The numbers are drawn from plausible ranges for a binder design campaign and
are correlated the way real metrics are: better interface energy tends to come
with more buried surface and higher predictor confidence, and a few designs are
deliberately bad so filters and the ranking exclusion have something to catch.

This is sample data for exercising the ranking and dashboard stages without a
multi-hour PyRosetta run. It is not experimental data and must not be used as
such.

Usage:
    python examples/generate_demo_metrics.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

N_DESIGNS = 60
SEED = 20260107
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def main() -> None:
    rng = np.random.default_rng(SEED)

    # One latent "quality" variable drives the correlations between metrics.
    quality = rng.beta(2.2, 2.0, N_DESIGNS)
    binder_len = rng.integers(65, 151, N_DESIGNS)

    dG = -8 - 34 * quality + rng.normal(0, 3.0, N_DESIGNS)
    dSASA = 500 + 1500 * quality + rng.normal(0, 140, N_DESIGNS)
    sc = np.clip(0.48 + 0.30 * quality + rng.normal(0, 0.035, N_DESIGNS), 0, 1)
    packstat = np.clip(0.44 + 0.30 * quality + rng.normal(0, 0.045, N_DESIGNS), 0, 1)
    n_hbonds = np.maximum(0, np.round(1 + 11 * quality + rng.normal(0, 1.6, N_DESIGNS)))
    unsat = np.maximum(0, np.round(7 - 6 * quality + rng.normal(0, 1.3, N_DESIGNS)))
    i_ptm = np.clip(0.33 + 0.58 * quality + rng.normal(0, 0.055, N_DESIGNS), 0, 1)
    ipsae = np.clip(0.22 + 0.62 * quality + rng.normal(0, 0.06, N_DESIGNS), 0, 1)
    plddt = np.clip(0.62 + 0.30 * quality + rng.normal(0, 0.04, N_DESIGNS), 0, 1)
    binder_plddt = np.clip(0.66 + 0.28 * quality + rng.normal(0, 0.05, N_DESIGNS), 0, 1)
    surface_hydro = np.clip(0.36 - 0.13 * quality + rng.normal(0, 0.035, N_DESIGNS), 0, 1)
    n_iface = np.maximum(4, np.round(9 + 16 * quality + rng.normal(0, 2.4, N_DESIGNS)))

    # Three relaxation artefacts with a positive interface energy. The ranking
    # exclusion in the config drops them before normalisation.
    dG[rng.choice(N_DESIGNS, 3, replace=False)] = rng.uniform(1.5, 9.0, 3)

    seeds = rng.integers(1000, 999999, N_DESIGNS)
    frame = pd.DataFrame(
        {
            "design": [
                f"demo_l{length}_s{seed}_mpnn{rng.integers(1, 21)}_model{rng.integers(1, 6)}"
                for length, seed in zip(binder_len, seeds, strict=False)
            ],
            "binder_len": binder_len,
            "clashes_pre_relax": np.maximum(
                0, np.round(6 - 5 * quality + rng.normal(0, 1.8, N_DESIGNS))
            ),
            "dG": dG,
            "dSASA": dSASA,
            "sc": sc,
            "packstat": packstat,
            "n_hbonds": n_hbonds,
            "unsat_hbonds": unsat,
            "n_interface_res": n_iface,
            "interface_hydro": np.clip(
                42 + 16 * quality + rng.normal(0, 7, N_DESIGNS), 0, 100
            ),
            "binder_sasa": 3600 + 42 * binder_len + rng.normal(0, 320, N_DESIGNS),
            "surface_hydro": surface_hydro,
            "binder_score": -1.9 * binder_len + rng.normal(0, 26, N_DESIGNS),
            "binder_helix_pct": np.clip(
                32 + 46 * quality + rng.normal(0, 13, N_DESIGNS), 0, 100
            ),
            "Average_i_pTM": i_ptm,
            "Average_ipSAE": ipsae,
            "Average_pLDDT": plddt,
            "Average_Binder_pLDDT": binder_plddt,
            "Sequence": [
                "".join(rng.choice(list(AMINO_ACIDS), size=int(length)))
                for length in binder_len
            ],
        }
    )

    frame["dG_dSASA_ratio"] = frame["dG"] / frame["dSASA"] * 100
    frame["hbond_pct"] = frame["n_hbonds"] / frame["n_interface_res"] * 100
    frame["unsat_hbonds_pct"] = frame["unsat_hbonds"] / frame["n_interface_res"] * 100
    frame["interface_sasa_pct"] = frame["dSASA"] / frame["binder_sasa"] * 100
    frame["design_base"] = frame["design"].str.replace(r"_model\d+$", "", regex=True)

    output = Path(__file__).parent / "demo_metrics.csv"
    frame.round(4).to_csv(output, index=False)
    print(f"wrote {output} ({len(frame)} synthetic designs)")


if __name__ == "__main__":
    main()
