# Contributing

Thanks for taking a look. Issues and pull requests are both welcome.

## Development setup

```bash
git clone https://github.com/jufercar/binderscope
cd binderscope
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
pytest
```

The test suite needs neither PyRosetta nor any experimental data, so it runs
anywhere in about a second. If you are changing the scoring stage you will need
a PyRosetta installation as well — see the README for the optional
dependencies.

## Checks

CI runs these on Python 3.10 through 3.13, and so should you before opening a
pull request:

```bash
ruff check .                              # lint and import order
pytest -q                                 # test suite
binderscope rank configs/demo.yaml        # the documented commands still work
binderscope dashboard configs/demo.yaml
```

`ruff check --fix .` handles import order and most style issues on its own.

## Reporting a problem

For anything hit on a real design campaign, please include the config file —
with confidential paths and target details removed — and what the stage
printed. Most failures come down to a chain id or a residue mapping that does
not match the structures, and the config plus the error message is usually
enough to see which.

## What goes where

Target-specific values never belong in code. Residue mappings, chain ids,
thresholds, weights and file paths are configuration, which is what lets a run
be reproduced from one YAML file — and lets the method be published while an
unpublished construct stays in a config that is not committed. If you find
yourself hard-coding a residue number, it probably wants a config key.

The modules are described in the README's project structure section.

## Adding a metric

1. Compute it where it belongs: `rosetta.py` for anything needing a relaxed
   pose, `spatial.py` for anything in the reference frame, `interface.py` for
   composition, `geometry.py` for pure coordinate maths.
2. Return it in the metric dictionary. Column names are part of the public
   interface, because configs refer to them — spatial metrics derive theirs
   from `clash_targets` labels, so new ones should follow that pattern.
3. Make it selectable from a config under `ranking.metrics` and `filters`
   rather than wiring it into the score.
4. Add it to `AXIS_HINTS` in `dashboard/builder.py` if it makes sense as a
   radar axis, with a short label and the direction that counts as better.
5. Document it in the README's metrics list.

## Adding a reversion strategy

Add it to `REVERSION_STRATEGIES` in `config.py` and handle it in
`reversion.py`. Please include a test that puts the reference structure in a
different coordinate frame from the design: a strategy that silently assumes
they share one produces side chains tens of ångström from their own backbone,
which is how the `donor` strategy was originally broken. `tests/conftest.py`
has a `reference_pdb_other_frame` fixture for exactly this.

## Tests

Tests build their own structures rather than shipping real ones — see
`tests/conftest.py` for the coarse fixtures used in superposition tests, and
`tests/test_sequences.py` for one with realistic peptide geometry. Prefer
asserting a property that would actually break (a bond length, a column name,
a count that must not saturate) over asserting whatever the code currently
returns.

## Commits

Conventional Commits (`fix:`, `feat:`, `docs:`, `test:`, `chore:`). Please say
in the body what breaks and why a change matters, not just what changed.

## Licence

Contributions are accepted under the MIT licence of this project.
