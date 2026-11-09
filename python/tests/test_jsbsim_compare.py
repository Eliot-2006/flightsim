"""JSBSim comparison harness checks (PRD 8.6, plan Task 6.3).

Pure-Python tests run everywhere (the CI `python` job has no C++ build); the end-to-end test needs build/release/flightsim.
"""
import csv
import shutil
from pathlib import Path

import numpy as np
import pytest

from jsbsim_compare import compare, gen_xml, run_jsbsim

REPO = Path(__file__).resolve().parents[2]
CLI = REPO / "build" / "release" / "flightsim"


def cols(t, **signals):
    return {"t_s": np.asarray(t, float), **{k: np.asarray(v, float) for k, v in signals.items()}}


def test_metric_is_difference_over_jsbsim_peak():
    t = np.linspace(0, 1, 101)
    jsb = cols(t, q_dps=10 * np.sin(2 * np.pi * t))
    ours = cols(t, q_dps=10 * np.sin(2 * np.pi * t) + 0.1 * t)  # worst difference 0.1 at t = 1
    assert compare.metric(ours, jsb, "q_dps") == pytest.approx(0.1 / 10, rel=1e-9)
    assert compare.metric(jsb, jsb, "q_dps") == 0.0


def test_metric_uses_deviation_from_each_own_trim():
    t = np.linspace(0, 1, 11)
    jsb = cols(t, h_ft=5000 + t**2)
    ours = cols(t, h_ft=7000 + t**2)  # different trim value, identical response
    assert compare.metric(ours, jsb, "h_ft") == 0.0


def test_metric_not_excited_signal_is_none():
    t = np.linspace(0, 1, 11)
    assert compare.metric(cols(t, p_dps=t * 1e-9), cols(t, p_dps=t * 1e-12), "p_dps") is None


def test_metric_subtracts_hold_drift_from_jsbsim():
    t = np.linspace(0, 1, 11)
    response, drift = np.sin(t), 0.05 * t**2
    ours, jsb, hold = cols(t, h_ft=response), cols(t, h_ft=response + drift), cols(t, h_ft=drift)
    assert compare.metric(ours, jsb, "h_ft") > 0.01
    assert compare.metric(ours, jsb, "h_ft", hold) == pytest.approx(0.0, abs=1e-12)


def test_variants_change_exactly_one_factor():
    base = (REPO / "scenarios/jsbsim/elev_doublet_h10000_dt2400.toml").read_text()
    g = compare.variant_text(base, "gravity_constant")
    assert "R = inf" in g and "R = 6378137" not in g and 'atmosphere = "us76"' in g
    a = compare.variant_text(base, "atmosphere_sl_adc")
    assert 'atmosphere = "sl_adc"' in a and "R = 6378137" in a
    jt = {"alpha_deg": 3.0, "theta_deg": 3.0, "throttle": 0.5, "de_deg": -0.6}
    t = compare.variant_text(base, "jsbsim_trim_state", jt)
    assert "from_trim = false" in t and "throttle = 0.5" in t and "de_deg = -0.6" in t
    x13 = [float(x) for x in t.split("x13 = [")[1].split("]")[0].split(",")]
    assert len(x13) == 13 and x13[12] == pytest.approx(64.94 * 0.5)  # power from the low-slope gearing segment
    assert x13[11] == pytest.approx(10000 * 0.3048) and x13[0] == pytest.approx(502 * 0.3048)


def test_our_planet_matches_jsbsim_planet_radius():
    import tomllib

    r = tomllib.loads((REPO / "scenarios/jsbsim/trim_h0.toml").read_text())["aircraft"]["gravity"]["R"]
    assert r == pytest.approx(gen_xml.R_FT * 0.3048, rel=1e-12)


def test_kinematic_curvature_estimate_matches_jsbsim_hold_at_5s():
    hold = compare.read_csv(compare.RESULTS / "trim_h10000_dt2400.csv")
    est = compare.analytic_curvature(502 * 0.3048, gen_xml.R_FT * 0.3048, 10000 * 0.3048, 5.0)
    i = 500  # 5 s at 0.01 s
    # 5 s is inside the kinematic regime (aerodynamic feedback is still < 1 % of the drift); 3 % allows the residual trim accel.
    assert est["h_ft"] == pytest.approx(hold["h_ft"][i] - hold["h_ft"][0], rel=0.03)
    assert est["theta_deg"] == pytest.approx(hold["theta_deg"][i] - hold["theta_deg"][0], rel=0.03)
    assert est["V_fps"] == pytest.approx(hold["V_fps"][i] - hold["V_fps"][0], rel=0.03)


def test_committed_reference_csvs_share_one_grid():
    for name, _ in compare.reference_runs():
        d = compare.read_csv(compare.RESULTS / f"{name}.csv")
        assert list(d) == run_jsbsim.COLS
        assert len(d["t_s"]) == 2001 and np.abs(np.diff(d["t_s"]) - 0.01).max() < 1e-9, name


def test_csv_matches_flags_a_perturbed_column(tmp_path):
    src = compare.RESULTS / "elev_doublet_h10000_dt2400.csv"
    shutil.copy(src, tmp_path / "same.csv")
    assert compare.csv_matches(tmp_path / "same.csv", src) == []
    rows = list(csv.reader(open(src)))
    rows[1000][rows[0].index("q_dps")] = repr(float(rows[1000][rows[0].index("q_dps")]) * (1 + 1e-6))
    csv.writer(open(tmp_path / "bad.csv", "w", newline="")).writerows(rows)
    assert [k for k, _ in compare.csv_matches(tmp_path / "bad.csv", src)] == ["q_dps"]


@pytest.mark.skipif(not CLI.exists(), reason="release CLI not built")
def test_end_to_end_criteria_and_budget_structure(tmp_path):
    res = compare.compute(CLI, tmp_path)
    assert all(r["pass"] for r in res["J1"]) and all(r["pass"] for r in res["J2"])
    assert all(r["pass"] for r in res["J3"].values()), {c: r["worst"] for c, r in res["J3"].items()}
    # the J1 w-dot offset is the density mismatch, to 1 % of itself
    for r in res["J1"]:
        assert r["jsbsim"]["wdot_fps2"] == pytest.approx(r["wdot_from_rho_fps2"], rel=0.01)
    assert set(res["budget"]) >= {"jsbsim_dt_halving", "our_dt_halving", "gravity_constant", "atmosphere", "curvature",
                                  "residual_trim_mismatch", "sampling_alignment", "table_edges", "all_explained"}
    # our own discretization is negligible next to the observed difference (PRD 8.6: "likewise negligible")
    ours = max(abs(v) for e in res["budget"]["our_dt_halving"].values() for v in e["self"].values() if v is not None)
    assert ours < 1e-6
    # ours column names match the JSBSim CSV (comparison convenience, plan Task 6.1)
    header = (tmp_path / "work/elev_doublet_ours/trajectory.csv").read_text().splitlines()[0].split(",")
    assert set(run_jsbsim.COLS) - {"throttle"} <= set(header) | {"throttle"}
