"""Report tooling tests. The fixture metrics.json is a copy of a real `out/metrics.json`; refresh it with
`cp ../out/metrics.json tests/fixtures/metrics.json` after the demo whenever the template gains a key."""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from flightsim_analysis import metrics, plots, report

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path(__file__).parent / "fixtures" / "metrics.json"
PLACEHOLDER = re.compile(r"\{\{.*?\}\}")


@pytest.mark.parametrize("template", ["results.template.md", "README.template.md"])
def test_render_fills_every_placeholder_in_real_template(tmp_path, template):
    m = json.loads(FIXTURE.read_text())
    out = tmp_path / "rendered.md"
    report.write_results(m, ROOT / "docs" / template, out)
    text = out.read_text()
    assert not PLACEHOLDER.search(text)
    assert "{{" not in text and "}}" not in text


def test_render_formats_values():
    assert report.render("p={{a.b:.2f}} n={{n}} s={{s}}", {"a.b": 1.23456, "n": 3, "s": "x"}) == "p=1.23 n=3 s=x"


def test_render_missing_key_raises():
    with pytest.raises(KeyError, match="nope"):
        report.render("{{nope}}", {})


def test_sign_changes_and_drift_helpers():
    assert metrics.sign_changes([1.0, 0.5, -0.2, -1.0, 0.3]) == 2
    assert metrics.max_rel_drift([2.0, 2.0, 2.2]) == pytest.approx(0.1)
    assert metrics.max_vec_rel_drift([[3, 4], [3, 4], [3, 4.5]]) == pytest.approx(0.5 / 5)


def test_figure_registry_has_required_figures():
    ids = {f.id for f in plots.FIGURES}
    assert {"convergence", "work-precision", "drift", "qnorm", "r4", "t4"} <= ids
    for f in plots.FIGURES:
        assert f.title and f.caption and f.test_ref  # every figure states what it shows and which test asserts it


CLI = ROOT / "build" / "release" / "flightsim"


@pytest.mark.skipif(not CLI.exists(), reason="release CLI not built")
def test_full_pipeline_builds_every_figure(tmp_path):
    from flightsim_analysis import run

    run.main(["--out", str(tmp_path / "out"), "--bench-min-time", "0.01s"])
    m = metrics.build(tmp_path / "out")
    html = report.build_html(tmp_path / "out", m)
    for f in plots.FIGURES:
        assert f'id="fig-{f.id}"' in html
    assert "<script src" not in html and "<link " not in html  # self-contained: Plotly JS inlined
    shutil.rmtree(tmp_path / "out")


# WCAG 1.4.11 non-text contrast: every series colour must reach 3:1 against the surface it is drawn on.
LIGHT_SURFACE, DARK_SURFACE = "#fcfcfb", "#1a1a19"  # keep equal to --bg in report.CSS
NON_TEXT_MIN = 3.0  # WCAG 2.1 SC 1.4.11


def _luminance(hex_color: str) -> float:
    c = [int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5)]
    lin = [x / 12.92 if x <= 0.03928 else ((x + 0.055) / 1.055) ** 2.4 for x in c]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def contrast(a: str, b: str) -> float:
    hi, lo = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


@pytest.mark.parametrize("colors,surface", [(plots.LIGHT, LIGHT_SURFACE), (report.DARK, DARK_SURFACE)])
def test_series_colours_meet_non_text_contrast(colors, surface):
    for c in colors + [plots.START_MARKER]:
        assert contrast(c, surface) >= NON_TEXT_MIN, f"{c} on {surface}: {contrast(c, surface):.2f}"


def test_surfaces_match_css():
    assert f"--bg:{LIGHT_SURFACE}" in report.CSS and f"--bg:{DARK_SURFACE}" in report.CSS
