"""J4 (PRD 8.6 appendix): JSBSim's stock f16 trimmed at the matched condition. Trim only; no pass/fail on agreement."""
import math
from pathlib import Path

import pytest

from jsbsim_compare import stock_f16

REPO = Path(__file__).resolve().parents[2]
CLI = REPO / "build" / "release" / "flightsim"


@pytest.fixture(scope="module", params=[0, 10000])
def trim(request):
    return stock_f16.stock_trim(request.param)


def test_trim_is_a_steady_state(trim):
    # JSBSim's own trim tolerance leaves accelerations of ~1e-4 (ft/s^2, rad/s^2); these bounds are 10x that.
    assert abs(trim["udot_fps2"]) < 1e-3 and abs(trim["wdot_fps2"]) < 1e-3 and abs(trim["qdot_rps2"]) < 1e-4
    assert trim["V_fps"] == pytest.approx(502.0, abs=1e-6)


def test_clean_configuration(trim):
    # The stock file starts gear-down, which adds a drag coefficient of ~1.8e3 (in its units); trim must be gear-up.
    assert trim["gear_pos_norm"] == 0.0 and trim["CDgear"] == 0.0
    assert trim["elevator_cmd_norm"] == 0.0  # stick centred: the pitch trim carries the trim


def test_elevator_mapping_is_the_surface_position(trim):
    # fcs/elevator-pos-rad is what the aero tables read; positive = trailing edge down, as in S&L.
    assert trim["elevator_deg"] == pytest.approx(math.degrees(trim["elevator_rad"]), rel=1e-12)
    assert abs(trim["elevator_deg"]) < 25.0  # +-0.436 rad range in the file
    assert trim["throttle_pos_norm"] == pytest.approx(2.0 * trim["throttle_cmd_norm"], rel=1e-12)


def test_stock_mass_properties_are_not_the_sl_values(trim):
    # 17,400 empty + 230 pilot + 3,000 fuel (f16.xml); S&L uses 20,500 lbf.
    assert trim["weight_lbf"] == pytest.approx(20630.0)


def test_altitude_changes_trim_the_right_way():
    lo, hi = stock_f16.stock_trim(0), stock_f16.stock_trim(10000)
    assert hi["alpha_deg"] > lo["alpha_deg"]  # lower dynamic pressure at the same V needs more lift coefficient
    assert hi["thrust_lbf"] < lo["thrust_lbf"] * 1.5  # sanity: same order of magnitude


@pytest.mark.skipif(not CLI.exists(), reason="release CLI not built")
def test_compute_pairs_ours_with_stock(tmp_path):
    rows = stock_f16.compute(CLI, tmp_path)
    assert [r["h_ft"] for r in rows] == [0, 10000]
    for r in rows:
        assert set(r["ours"]) >= {"alpha_deg", "elevator_deg", "throttle", "thrust_lbf"}
        assert r["difference"]["alpha_deg"] == pytest.approx(r["stock"]["alpha_deg"] - r["ours"]["alpha_deg"])
