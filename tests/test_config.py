from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from binderscope.config import ConfigError, load_config

MINIMAL = {
    "name": "unit",
    "designs": "./designs",
    "output": "./analysis",
    "chains": {"target": "A", "binder": "B"},
}

CONFIG_DIR = Path(__file__).resolve().parents[1] / "configs"


def write(tmp_path: Path, data: dict) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


def test_minimal_config_loads(tmp_path):
    config = load_config(write(tmp_path, MINIMAL))
    assert config.name == "unit"
    assert config.target_chain == "A"
    assert config.reference is None
    assert not config.needs_reversion


def test_relative_paths_resolve_against_the_config_file(tmp_path):
    config = load_config(write(tmp_path, MINIMAL))
    assert config.designs_dir == (tmp_path / "designs").resolve()
    assert config.output_dir == (tmp_path / "analysis").resolve()


def test_derived_output_paths_are_namespaced_by_run(tmp_path):
    config = load_config(write(tmp_path, MINIMAL))
    assert config.metrics_csv.name == "unit_metrics.csv"
    assert config.ranked_csv.name == "unit_ranked.csv"
    assert config.dashboard_html.name == "unit_dashboard.html"


def test_missing_required_key_is_reported(tmp_path):
    data = {k: v for k, v in MINIMAL.items() if k != "designs"}
    with pytest.raises(ConfigError, match="missing required key: 'designs'"):
        load_config(write(tmp_path, data))


def test_absent_file_is_reported():
    with pytest.raises(ConfigError, match="config file not found"):
        load_config("/nonexistent/config.yaml")


def test_reversion_requires_a_native_target(tmp_path):
    data = dict(MINIMAL)
    data["reference"] = {
        "structure": "./ref.pdb",
        "target_chain": "D",
        "segments": [{"input": [1, 10], "reference_start": 101}],
        "revert_mutations": [
            {"position": 5, "wild_type": "ALA", "reference_position": 105}
        ],
    }
    with pytest.raises(ConfigError, match="requires reference.native_target"):
        load_config(write(tmp_path, data))


def test_invalid_amino_acid_code_is_rejected(tmp_path):
    data = dict(MINIMAL)
    data["reference"] = {
        "structure": "./ref.pdb",
        "target_chain": "D",
        "native_target": "./native.pdb",
        "revert_mutations": [
            {"position": 5, "wild_type": "XYZ", "reference_position": 105}
        ],
    }
    with pytest.raises(ConfigError, match="not a three-letter amino acid code"):
        load_config(write(tmp_path, data))


def test_unsupported_filter_operator_is_rejected(tmp_path):
    data = dict(MINIMAL)
    data["filters"] = [{"col": "dG", "op": "<=", "value": -10}]
    with pytest.raises(ConfigError, match="unsupported operator"):
        load_config(write(tmp_path, data))


def test_interface_is_derived_from_the_chains(tmp_path):
    data = dict(MINIMAL)
    data["chains"] = {"target": "C", "binder": "D"}
    config = load_config(write(tmp_path, data))
    assert config.rosetta.interface == "C_D"


def test_interface_disagreeing_with_the_chains_is_rejected(tmp_path):
    # Would otherwise have InterfaceAnalyzerMover score chains that are absent.
    data = dict(MINIMAL)
    data["chains"] = {"target": "C", "binder": "D"}
    data["rosetta"] = {"interface": "A_B"}
    with pytest.raises(ConfigError, match="rosetta.interface is 'A_B'"):
        load_config(write(tmp_path, data))


def test_interface_agreeing_with_the_chains_is_accepted(tmp_path):
    data = dict(MINIMAL)
    data["rosetta"] = {"interface": "A_B"}
    assert load_config(write(tmp_path, data)).rosetta.interface == "A_B"


def test_reversion_strategy_defaults_to_truncate(tmp_path):
    assert load_config(write(tmp_path, MINIMAL)).reversion_strategy == "truncate"


def test_unknown_reversion_strategy_is_rejected(tmp_path):
    data = dict(MINIMAL)
    data["reversion_strategy"] = "magic"
    with pytest.raises(ConfigError, match="reversion_strategy must be one of"):
        load_config(write(tmp_path, data))


def test_reversion_strategy_under_pymol_is_rejected(tmp_path):
    # It used to be read from there, which made no sense and was easy to miss.
    data = dict(MINIMAL)
    data["pymol"] = {"reversion_strategy": "donor"}
    with pytest.raises(ConfigError, match="belongs at the top level"):
        load_config(write(tmp_path, data))


def test_reference_requires_structure_and_chain(tmp_path):
    data = dict(MINIMAL)
    data["reference"] = {"target_chain": "D"}
    with pytest.raises(ConfigError, match="requires 'structure' and 'target_chain'"):
        load_config(write(tmp_path, data))


@pytest.mark.parametrize(
    "filename", ["demo.yaml", "pdl1_example.yaml", "engineered_target_template.yaml"]
)
def test_shipped_configs_are_valid(filename):
    config = load_config(CONFIG_DIR / filename)
    assert config.name
    assert config.rank_metrics


def test_geometry_options_reach_the_scorer():
    # The scorer used to be built with RosettaOptions alone, so every
    # documented geometry threshold was silently ignored.
    import inspect

    from binderscope.rosetta import InterfaceScorer, build_scorer

    assert "geometry" in inspect.signature(InterfaceScorer.__init__).parameters
    assert "geometry=config.geometry" in inspect.getsource(build_scorer)

    source = inspect.getsource(InterfaceScorer.score)
    assert "self.geometry.prerelax_clash_cutoff" in source
    assert "self.geometry.interface_cutoff" in source


def test_engineered_template_exercises_the_full_schema():
    config = load_config(CONFIG_DIR / "engineered_target_template.yaml")
    assert config.needs_reversion
    assert config.reference is not None
    assert len(config.reference.segments) == 3
    assert len(config.reference.revert_mutations) == 4
    assert {t.label for t in config.reference.clash_targets} == {
        "target",
        "native_partner",
    }
    assert config.reference.linker_positions == {74, 75, 76, 84, 85}
