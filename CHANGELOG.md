# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — 2026-10-07

First release. A post-processing pipeline for protein binder design campaigns,
refactored from a single-target analysis notebook into a configuration-driven
tool, so that the method can be reproduced and published independently of the
data it was developed on.

### Added

- **Wild-type reversion** (`binderscope revert`). Restores engineered
  substitutions in a design target before any metric is computed, so scores
  describe the biologically relevant interface rather than the engineered
  surface. Two strategies: `truncate`, which keeps the backbone and `CB` and
  lets Rosetta rebuild the side chain, and `donor`, which transplants the
  reference rotamer.
- **PyRosetta interface scoring** (`binderscope score`). Chi-only `FastRelax`,
  `InterfaceAnalyzerMover`, `BuriedUnsatHbonds` via DAlphaBall, total energy,
  SASA and layer metrics, interface composition and secondary structure.
- **Placement geometry.** Superposes a complex back into a reference frame
  through a configurable residue map and counts binder contacts against chosen
  reference chains, with voxelised volume overlap — answering where the binder
  ended up, which interface energy cannot.
- **Predictor merge** (`binderscope merge`). Joins the AlphaFold2 statistics
  BindCraft already writes, tolerating version differences in the columns.
- **Multi-criteria ranking** (`binderscope rank`). Weighted, normalised scoring
  with per-metric contributions reported alongside the total, and exclusions
  applied before normalisation so a single artefact cannot stretch the scale.
- **Self-contained dashboard** (`binderscope dashboard`). One HTML file with
  the data inlined: context panel generated from the config, live filters, a
  comparison radar, FASTA export and a ranking view.
- **Sequence export** (`binderscope sequences`) to FASTA and CSV.
- **Resumable stages.** Results are appended as they are produced and a restart
  skips completed work; a design that fails is recorded and the run continues.
- Optional per-design PyMOL sessions (`revert --sessions`).
- Three annotated configurations: a synthetic demo needing no dependencies, the
  public PD-L1 benchmark, and a template for engineered targets.
- A walkthrough notebook, synthetic sample data with its generator, and a test
  suite that requires neither PyRosetta nor experimental data.

### Notes

The PyRosetta stage cannot be covered by the test suite, so it is verified by
running it: `score` has been exercised on real two-chain complexes with
DAlphaBall and DSSP active, and the placement-geometry output reproduces the
reference implementation this was refactored from to all printed decimals.

Several defects were found and fixed during that verification, before this
release: side-chain grafting ignored the reference's own coordinate frame and
then its per-residue backbone, distorting bonds that chi-only relax cannot
repair; `BuriedUnsatHbonds` never loaded because of a wrong call signature
hidden by a silent fallback; and the documented geometry thresholds never
reached the scorer. None of these ever shipped.

[Unreleased]: https://github.com/jufercar/binderscope/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/jufercar/binderscope/releases/tag/v0.1.0
