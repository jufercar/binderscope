# Example data

## `demo_metrics.csv` — synthetic

Sixty designs of **randomly generated** metrics, produced by
`generate_demo_metrics.py`. A single latent quality variable drives the
correlations, so better interface energy tends to come with more buried surface
and higher predictor confidence, the way real metrics behave. Three designs are
given a positive interface energy on purpose, so the ranking exclusion has
something to catch.

**This is not experimental data.** It exists so that

```bash
binderscope rank      configs/demo.yaml --contributions
binderscope dashboard configs/demo.yaml
```

run in seconds without PyRosetta, a GPU or a design campaign. Do not use it for
anything else. Regenerate it with:

```bash
python examples/generate_demo_metrics.py
```

The seed is fixed, so the file is reproducible.

## Real structures

None are included. Design PDBs, trajectories and model weights are large and
often unpublished, so they stay out of the repository — point a config at your
own run directory instead. To reproduce a real run end to end, run the
BindCraft PD-L1 example and use `configs/pdl1_example.yaml`.
