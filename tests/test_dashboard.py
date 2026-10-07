from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest

from binderscope.config import load_config
from binderscope.dashboard import build_dashboard
from binderscope.merge import merge_predictor_stats, split_design_name

ROOT = Path(__file__).resolve().parents[1]
DEMO_CONFIG = ROOT / "configs" / "demo.yaml"
DEMO_METRICS = ROOT / "examples" / "demo_metrics.csv"


@pytest.fixture
def demo_frame() -> pd.DataFrame:
    return pd.read_csv(DEMO_METRICS)


def build(tmp_path, frame, config_path=DEMO_CONFIG):
    config = load_config(config_path)
    return build_dashboard(frame, config, tmp_path / "dashboard.html")


def test_dashboard_is_written(tmp_path, demo_frame):
    output = build(tmp_path, demo_frame)
    assert output.is_file()
    assert output.stat().st_size > 20_000


def test_every_placeholder_is_filled(tmp_path, demo_frame):
    html = build(tmp_path, demo_frame).read_text()
    assert not re.findall(r"__[A-Z_]+__", html)


def test_data_is_inlined_so_the_file_stands_alone(tmp_path, demo_frame):
    html = build(tmp_path, demo_frame).read_text()
    match = re.search(r"const DATA\s*=\s*(\[.*?\]);\n", html, re.S)
    assert match
    records = json.loads(match.group(1))
    assert len(records) == len(demo_frame)
    assert "design" in records[0]


def test_context_panels_come_from_the_configuration(tmp_path, demo_frame):
    html = build(tmp_path, demo_frame).read_text()
    assert "Run: demo" in html
    assert "Multi-criteria ranking" in html
    # No reference section in the demo config, so these panels must be absent.
    assert "Target mapping" not in html
    assert "Reverted substitutions" not in html


def test_reference_configuration_adds_its_panels(tmp_path, demo_frame):
    html = build(
        tmp_path, demo_frame, ROOT / "configs" / "engineered_target_template.yaml"
    ).read_text()
    assert "Target mapping" in html
    assert "Reverted substitutions" in html
    assert "136</b> mapped residues" in html


def test_ranking_weights_reach_the_page(tmp_path, demo_frame):
    html = build(tmp_path, demo_frame).read_text()
    match = re.search(r"const RANK_METRICS=(\[.*?\]);", html, re.S)
    assert match
    metrics = json.loads(match.group(1))
    assert {"col", "label", "higher", "w"} <= set(metrics[0])
    assert any(m["col"] == "dG" and m["higher"] is False for m in metrics)


def test_empty_table_is_refused(tmp_path):
    with pytest.raises(ValueError, match="empty metrics table"):
        build(tmp_path, pd.DataFrame())


TEMPLATE = ROOT / "src" / "binderscope" / "dashboard" / "template.html"


def test_template_is_target_agnostic():
    # Guards the property that makes a run publishable: the shipped template
    # describes no particular target. Anything target-specific must arrive
    # through a placeholder at build time, never be baked into the HTML.
    template = TEMPLATE.read_text()
    for placeholder in ("__TITLE__", "__CONTEXT_HTML__", "__DATA__", "__RANK_METRICS__"):
        assert placeholder in template

    # No residue-range or chain-selection literals, which is how construct
    # details leaked into the notebook this template came from.
    assert not re.search(r"\bresi\s+\d+\s*-\s*\d+", template)
    assert not re.search(r"\bchain\s+[A-Z]\b", template)


def test_template_carries_no_embedded_run_data():
    template = TEMPLATE.read_text()
    assert "const DATA   = __DATA__;" in template
    # A design name would mean a previous run's data was committed with it.
    assert not re.search(r"_mpnn\d+_model\d+", template)


# ── merge helpers ────────────────────────────────────────────────────────────
@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("design_l50_s123_mpnn1_model2", ("design_l50_s123_mpnn1", 2)),
        ("design_l50_s123_mpnn1", ("design_l50_s123_mpnn1", None)),
        ("plain", ("plain", None)),
    ],
)
def test_model_suffix_is_split_off(name, expected):
    assert split_design_name(name) == expected


def test_merge_is_a_no_op_without_statistics(demo_frame):
    assert merge_predictor_stats(demo_frame, None) is demo_frame


def test_merge_joins_on_the_design_base(tmp_path):
    metrics = pd.DataFrame(
        {"design": ["d_model1", "d_model2"], "design_base": ["d", "d"], "dG": [-20.0, -25.0]}
    )
    stats = tmp_path / "stats.csv"
    pd.DataFrame({"Design": ["d"], "Average_i_pTM": [0.81], "Sequence": ["AAA"]}).to_csv(
        stats, index=False
    )
    merged = merge_predictor_stats(metrics, stats)
    assert list(merged["Average_i_pTM"]) == [0.81, 0.81]
    assert len(merged) == 2


def test_merge_reports_a_missing_statistics_file(demo_frame, tmp_path):
    with pytest.raises(FileNotFoundError, match="statistics file not found"):
        merge_predictor_stats(demo_frame, tmp_path / "nope.csv")
