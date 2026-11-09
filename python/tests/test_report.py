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
    assert {"trim.table", "linear.modes.table", "linear.matrices.table", "linear.ln2.ratio", "linear.ln3.table"} <= m.keys()
    assert {"ln3", "ln2"} <= {f.id for f in plots.FIGURES}
    assert m["linear.ln2.status"] == "pass"
    # Five Markdown columns plus the outer delimiters: embedded |beta| must be escaped.
    for row in m["linear.ln2.changes.table"].splitlines():
        assert len(re.findall(r"(?<!\\)\|", row)) == 6
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


@pytest.mark.skipif(not CLI.exists(), reason="release CLI not built")
@pytest.mark.parametrize("cg", ["030", "035"])
def test_aircraft_cli_json(tmp_path, cg):
    for command in ("trim", "linearize"):
        dest = tmp_path / f"{command}.json"
        subprocess.run([CLI, command, ROOT / f"scenarios/f16_trim_502_xcg{cg}.toml", "--out", dest], check=True)
        result = json.loads(dest.read_text())
        assert result["trim"]["ok"]
        # Task 5.1 convergence threshold, reported by C++.
        assert result["trim"]["residual"] < 1e-10
        assert len(result["trim"]["x17"]) == 17
        if command == "linearize":
            assert len(result["A"]) == 13 and all(len(row) == 13 for row in result["A"])
            assert len(result["B"]) == 13 and all(len(row) == 4 for row in result["B"])
            assert result["input_units"] == ["1", "rad", "rad", "rad"]
            if cg == "030":
                assert {m["name"] for m in result["modes"]} == {"phugoid", "short_period", "dutch_roll", "roll", "spiral"}
            else:
                assert any(m["real"] > 0 for m in result["modes"])
            assert len(result["sweep"]) == 9
            assert len(result["responses"]) == 2
            assert result["response_operating_point"]["h_m"] == 30.48
            if cg == "030":
                for response in result["responses"]:
                    assert response["max_abs_aileron"] == response["max_abs_rudder"] == 0
                    assert not response["out_of_domain"] and not response["cells_changed"]
                    assert len(response["cell_changes"]) == 6
                    assert all(change["ignored"] for change in response["cell_changes"])


def test_checks_parse_the_numbers_the_tests_print(tmp_path):
    # Fixtures are the real stdout of the A1/B1/B2/S1/W1/TR2/aircraft-convergence tests (see run.CHECKS).
    (tmp_path / "checks").mkdir()
    for name in ("validation", "unit"):
        shutil.copy(Path(__file__).parent / "fixtures" / f"checks_{name}.txt", tmp_path / "checks" / f"{name}.txt")
    c = metrics.checks(tmp_path)
    assert c["checks.a1.max_rel"] == pytest.approx(1.99905157675e-4)  # p_dot, rounded S&L inertia coefficients
    assert c["checks.a1.max_rel_excluding_rotational"] < 1e-4  # PRD 8.5 investigation threshold
    assert 3.9 < c["checks.b1.order"] < 4.1 and c["checks.b1.residual"] < 1e-6
    assert 2.0 < c["checks.aircraft.position_order"] < 2.5  # reduced order from piecewise-linear tables (PRD 2)
    assert c["checks.s1.symmetrized"] < 1e-12 and c["checks.w1.air_relative"] < 1e-9


def test_checks_fail_loudly_when_a_line_is_missing(tmp_path):
    (tmp_path / "checks").mkdir()
    (tmp_path / "checks" / "validation.txt").write_text("B1 dt=0.01 residual=1e-10\n")
    (tmp_path / "checks" / "unit.txt").write_text("")
    with pytest.raises(ValueError):
        metrics.checks(tmp_path)


def test_bench_table_puts_target_beside_measurement_and_states_a_miss():
    m = {key: limit / 2 for key, _, limit, *_ in metrics.BENCH_TARGETS}
    assert metrics.bench_f16_table(m)["bench.f16.misses"] == 0
    m["bench.f16.BM_F16Deriv.ns"] = 401.0  # PRD 11 target is 400 ns
    t = metrics.bench_f16_table(m)
    assert t["bench.f16.misses"] == 1 and "**missed**" in t["bench.f16.table"] and "≤ 400 ns" in t["bench.f16.table"]
