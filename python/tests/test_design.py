"""C5 is a loop-at-a-time requirement, including the other closed loops."""
import json
import tomllib

import numpy as np

from control_design import design


def _assert_close(a, b, path="gains"):
    assert type(a) is type(b) or (isinstance(a, (int, float)) and isinstance(b, (int, float))), path
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            _assert_close(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _assert_close(x, y, f"{path}[{i}]")
    elif isinstance(a, (int, float)) and not isinstance(a, bool):
        assert np.isclose(a, b, rtol=1e-9, atol=0.0), f"{path}: {a} != {b}"
    else:
        assert a == b, path


def test_committed_gains_meet_c5():
    gains = tomllib.loads(design.GAINS.read_text())
    plant = design.load_plant()
    report = design.evaluate(plant, gains)
    for loop in report.values():
        # PRD C5: 6 dB gain margin, 45 degrees phase margin.
        assert loop["gm_db"] >= 6
        assert loop["pm_deg"] >= 45
        # C5 also requires a factor-of-two reserve toward lower-gain boundaries.
        if loop["gain_reduction_margin_db"] is not None:
            assert loop["gain_reduction_margin_db"] >= 6
        assert loop["crossover_rad_s"] > 0
    # A return ratio's margins alone cannot certify an unstable aircraft.
    assert design.closed_loop_stable(plant, gains)
    # Required minimum nested-loop crossover separation (orchestrator: 3–5x).
    assert report["pitch_rate"]["crossover_rad_s"] >= 3 * report["pitch"]["crossover_rad_s"]
    assert report["pitch"]["crossover_rad_s"] >= 3 * report["altitude"]["crossover_rad_s"]


def test_reproducible_outputs(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    design.generate(first)
    design.generate(second)
    # Committed vs regenerated: numerically equal. Bitwise identity holds on one platform (asserted below), but
    # LAPACK/libm differ across platforms in the last digits, so the committed file is compared at 1e-9 relative.
    _assert_close(tomllib.loads((first / "autopilot_gains.toml").read_text()), tomllib.loads(design.GAINS.read_text()))
    for name in ("autopilot_gains.toml", "linear_steps.json", "margins.json"):
        assert (first / name).read_bytes() == (second / name).read_bytes()
    responses = json.loads((first / "linear_steps.json").read_text())
    for case, command in (("C1", 30.48), ("C3", 15.24)):
        response = responses[case]
        assert response["step_si"] == command  # Exact ft→m conversion at the reporting boundary.
        assert len(response["time_s"]) == len(response["altitude_m"])
        assert np.isfinite(response["airspeed_mps"]).all()


def test_linear_c1_headroom():
    gains = tomllib.loads(design.GAINS.read_text())
    response = design.step_predictions(design.load_plant(), gains)["C1"]
    time = np.array(response["time_s"])
    altitude = np.array(response["altitude_m"])
    speed = np.array(response["airspeed_mps"])
    step = response["step_si"]
    # Orchestrator headroom: 8% overshoot, 20 s settling, 4 ft/s speed excursion.
    assert max(altitude) <= 1.08 * step
    assert np.all(np.abs(altitude[time >= 20] - step) <= .02 * step)
    assert max(np.abs(speed)) <= 4 * .3048


def test_linear_c3_headroom():
    gains = tomllib.loads(design.GAINS.read_text())
    response = design.step_predictions(design.load_plant(), gains)["C3"]
    # Orchestrator headroom: 8% overshoot, 15 ft altitude excursion.
    assert max(response["airspeed_mps"]) <= 1.08 * response["step_si"]
    assert max(np.abs(response["altitude_m"])) <= 15 * .3048
