# binderscope

Structure-aware triage of protein binder designs.

A design campaign gives you hundreds of candidates and a table of numbers.
Choosing which ones to order is where the work actually is. `binderscope` is the
post-processing stage: it scores interfaces with PyRosetta, checks whether a
binder sits where it was meant to sit, ranks designs on criteria you declare,
and renders the result as a dashboard you can open, filter and hand to a
collaborator.

It reads [BindCraft](https://github.com/martinpacesa/BindCraft) output and works
with any pipeline that writes a two-chain complex per design.

```bash
binderscope rank      configs/demo.yaml --contributions
binderscope dashboard configs/demo.yaml
```

Those two commands run on synthetic sample data, need no PyRosetta and no GPU,
and produce a working dashboard in a couple of seconds.

## What it does that a metrics table does not

**Restores the target to wild type before scoring.** Design campaigns often run
against a modified target: surface residues swapped to steer the binder away
from an unwanted face, fragments of a structure spliced together with linkers.
Scoring against that construct measures an interface that does not exist in
nature. `binderscope` superposes the as-designed target onto each prediction,
restores the wild-type residues in place, and leaves the binder untouched — so
every metric describes the interface you actually care about.

**Asks where the binder ended up.** Interface energy says how good a complex is
on its own terms. It says nothing about placement. Superposing the complex back
into a reference structure and counting binder atoms that contact chosen
reference chains answers a different question: a binder penetrating its own
target is wrong however well it scores, and one occupying the same space as a
known natural partner is evidence it engages the intended surface. A voxelised
overlap volume quantifies that mimicry.

**Ranks on stated priorities, and shows its work.** No single metric picks a
binder; interface energy, shape complementarity, predictor confidence and
placement geometry each catch what the others miss. Every configured metric is
normalised in its own good direction, weighted and summed — and the per-metric
contributions are reported next to the total, because *why* a design ranks where
it does matters more than the rank. Weights live in the config file where a
reader can see and argue with them.

**Keeps the target out of the code.** Residue mappings, linker positions,
substitutions, reference chains, thresholds and weights are all configuration.
The package carries no knowledge of any particular target, which is what makes
a run reproducible by someone else — and what lets you publish the method while
an unpublished construct stays in a config file you do not commit.

## Install

```bash
git clone https://github.com/jcarrillo/binderscope
cd binderscope
pip install -e .
```

Python ≥ 3.10, plus Biopython, NumPy, pandas, SciPy and PyYAML.

Three things are optional, and only the stages that need them are affected:

| Component | Needed for | Without it |
|---|---|---|
| [PyRosetta](https://www.pyrosetta.org/) | `score` | `rank` and `dashboard` still run on an existing metrics CSV |
| `DAlphaBall.gcc` | buried unsatisfied H-bonds | falls back to InterfaceAnalyzer's own estimate |
| `dssp` | secondary structure | those columns are left empty |
| PyMOL | `revert --sessions` | sessions are skipped |

`DAlphaBall.gcc` and `dssp` ship with BindCraft, in its `functions/` directory.
PyRosetta is free for academic use but requires a licence; it is not bundled.

## Use

A run is one YAML file. Three are included:

- **`configs/demo.yaml`** — synthetic data, no dependencies, for seeing the output.
- **`configs/pdl1_example.yaml`** — the public PD-L1 benchmark target. The simplest
  real case: the target is the native protein, so no reversion is needed.
- **`configs/engineered_target_template.yaml`** — an annotated template for a
  chimeric target mapped onto a crystallographic reference, exercising the full
  schema. All values are placeholders.

```bash
binderscope revert    config.yaml [--sessions]   # rebuild the wild-type target
binderscope score     config.yaml                # PyRosetta + placement geometry
binderscope merge     config.yaml                # join BindCraft's predictor stats
binderscope rank      config.yaml --contributions
binderscope dashboard config.yaml
binderscope run       config.yaml                # all of the above
```

Stages are separate commands because their costs differ by orders of magnitude.
Scoring is hours of FastRelax; ranking and the dashboard are seconds and get
re-run constantly while thresholds are tuned. `score` and `revert` append
results as they go and skip completed work on restart, so an interrupted run
resumes rather than starting over, and one malformed PDB is recorded as a
failure instead of ending the run.

### Metrics

Per design, from PyRosetta (`InterfaceAnalyzerMover` after a chi-only
`FastRelax`): `dG`, `dSASA`, `dG_dSASA_ratio`, `sc`, `packstat`, `n_hbonds`,
`unsat_hbonds`, `n_interface_res`, `interface_hydro`, `binder_score`,
`binder_sasa`, `interface_sasa_pct`, `surface_hydro`, and helix/sheet/loop
percentages for the binder and its interface.

With a `reference` section: `rmsd_reference`, `n_superposition_anchors`, and per
reference chain `clash_<label>_hard`, `clash_<label>_soft`, plus
`V_binder_A3`, `V_inter_<label>_A3`, `V_inter_<label>_pct` where an overlap
volume was requested.

With `bindcraft_stats`: the AlphaFold2 columns BindCraft already wrote
(`Average_i_pTM`, `Average_ipSAE`, `Average_pLDDT`, MPNN scores, sequences and
the rest), joined per design.

### Dashboard

One self-contained HTML file with the data inlined — no server, no install, and
it survives being emailed or dropped in a supplement. Plotly is the only
external asset. Three tabs: a **Context** panel generated from the config
(target mapping, reverted substitutions, pipeline steps and ranking weights, so
the page documents its own provenance), an **Explorer** with live filters, a
comparison radar and FASTA export, and a **Ranking** view with per-metric
contribution breakdowns.

## Notes on method

**Reversion strategy.** The default, `truncate`, keeps the backbone and `CB`,
relabels the residue, and lets Rosetta rebuild the side chain with ideal
geometry before the chi-only relax repacks it — chemically consistent for any
substitution, including ones where the wild-type residue is larger than the
engineered one. The alternative, `donor`, transplants heavy-atom coordinates
from the reference structure, preserving the experimentally observed rotamer at
the cost of requiring that residue to be present and complete.

**Why chi-only relax.** The backbone comes from a structure predictor and is the
hypothesis under test. Letting FastRelax move it would improve the Rosetta score
while destroying the thing being measured. Side chains are repacked because
predicted rotamers are not reliable enough to score directly.

**Why exclusions happen before normalisation.** A single relaxation artefact
with a positive interface energy will stretch a min-max scale and compress every
real design into a narrow band. The `ranking.exclude` rules drop such rows
first; they are reported, not silently discarded.

## Sample data

`examples/demo_metrics.csv` is **synthetic**, generated by
`examples/generate_demo_metrics.py` with realistic correlations between metrics.
It exists so the ranking and dashboard stages can be exercised without a
multi-hour run. It is not experimental data and must not be used as such.

No structures are included: PDBs, trajectories and model weights are large and
frequently unpublished. Point a config at your own run instead — `.gitignore` is
set up to keep structures and run outputs out of the repository.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

66 tests covering residue mapping, clash and volume geometry, reversion against
synthetic structures, ranking and filter semantics, config validation and
dashboard rendering. None require PyRosetta or any experimental data.

## Relationship to BindCraft

This is a separate post-processing tool, not a fork. It consumes BindCraft
output and does not modify or redistribute it. If you use it, cite BindCraft
as well:

> Pacesa, M., Nickel, L., Schellhaas, C. *et al.* One-shot design of functional
> protein binders with BindCraft. *Nature* **646**, 483–492 (2025).
> https://doi.org/10.1038/s41586-025-09429-6

## Citing

See [`CITATION.cff`](CITATION.cff).

## Licence

MIT — see [`LICENSE`](LICENSE).
