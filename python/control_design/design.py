"""Reproduce M7 gains, loop-at-a-time margins, and C1/C3 linear predictions.

Run from python/: uv run python -m control_design.design
Requires a release CLI build. No nonlinear simulation is used for synthesis.
"""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import tomllib

import control as ct
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GAINS = ROOT / "scenarios/autopilot_gains.toml"
LOOPS = ("pitch_rate", "pitch", "altitude", "airspeed")
RATE_HZ = 50.0
# Loop-shaping choices, rad/s: 5x then 5x nominal cascade separation.
CROSSOVERS = (7.5, 1.5, 0.3, 0.12)
ZEROS = (0.5, 0.0, 0.01, 0.01)
# Resolve real crossings without ill-conditioned high-order polynomial conversion.
OMEGA = np.geomspace(1e-5, 1e4, 6000)


def load_plant():
    """CLI ideal-surface Jacobian + elevator lag + half-period Padé delays."""
    scenario = ROOT / "scenarios/f16_ap_trim.toml"
    text = scenario.read_text()
    settings = tomllib.loads(text)
    data_path = ROOT / settings["aircraft"]["data"]
    data = tomllib.loads(data_path.read_text())
    # The CLI explicitly rejects actuator states: preserve the approved scenario.
    ideal = text.replace("actuators = true", "actuators = false").replace(
        f'"{settings["aircraft"]["data"]}"', json.dumps(str(data_path))
    )
    with tempfile.TemporaryDirectory() as temp:
        source, dest = Path(temp) / "trim.toml", Path(temp) / "linear.json"
        source.write_text(ideal)
        subprocess.run(
            [str(ROOT / "build/release/flightsim"), "linearize", str(source), "--out", str(dest)],
            check=True, capture_output=True, text=True, cwd=ROOT,
        )
        model = json.loads(dest.read_text())
    if not model["trim"]["ok"] or model["trim"]["lookup_out_of_domain"]:
        raise ValueError("control trim must converge inside the data envelope")
    # [V, alpha, theta, q, h, P]; retain engine lag, omit lateral/navigation modes.
    ix = [0, 1, 4, 7, 11, 12]
    a = np.array(model["A"])[np.ix_(ix, ix)]
    b = np.array(model["B"])[np.ix_(ix, [1, 0])]
    c = np.eye(6)[[3, 2, 4, 0]]  # q, theta, h, V
    aircraft = ct.ss(a, b, c, np.zeros((4, 2)))
    s = ct.tf("s")
    delay = ct.tf(*ct.pade(1 / (2 * RATE_HZ), 1))
    lag = data["actuators"]["elevator"]["lag_rate"]
    physical = aircraft * ct.append(ct.ss(lag / (s + lag) * delay), ct.ss(delay))
    # u=[elevator,q_cmd,theta_cmd,throttle]; y=[q-q_cmd,theta-theta_cmd,h,V].
    # Negative diagonal feedback of C_i therefore gives the prescribed cascade.
    inputs = np.array([[1, 0, 0, 0], [0, 0, 0, 1]])
    direct = np.zeros((4, 4))
    direct[0, 1] = direct[1, 2] = -1
    return ct.ss(physical.A, physical.B @ inputs, physical.C, physical.D @ inputs + direct)


def controllers(gains):
    """Diagonal P/PI blocks; no unused integrator state in the pitch P loop."""
    blocks = []
    for name in LOOPS:
        g = gains[name]
        if g["ki"]:
            blocks.append(ct.ss(ct.tf([g["kp"], g["ki"]], [1, 0])))
        else:
            blocks.append(ct.ss([], [], [], [[g["kp"]]]))
    return ct.append(*blocks)


def return_ratio(plant, gains, index):
    """Break exactly one error signal; every other controller remains closed."""
    mask = np.eye(4)
    mask[index, index] = 0
    return ct.feedback(plant * controllers(gains), mask)[index, index]


def closed_loop(plant, gains):
    return ct.feedback(plant * controllers(gains), np.eye(4))


def closed_loop_stable(plant, gains):
    return bool(np.all(ct.poles(closed_loop(plant, gains)).real < 0))


def evaluate(plant, gains):
    report = {}
    for i, name in enumerate(LOOPS):
        ratio = return_ratio(plant, gains, i)
        gm, pm, _, wpc, wgc, _ = ct.stability_margins(ct.frd(ratio, OMEGA), returnall=True)
        boundaries = 20 * np.log10(gm)
        increasing = boundaries[boundaries > 0]
        decreasing = -boundaries[boundaries < 0]
        if not len(pm):
            raise ValueError(f"{name}: no gain crossover")
        # Upper and lower gain boundaries are distinct for nested feedback.
        report[name] = {
            "crossover_rad_s": float(wgc[np.argmin(pm)]),
            "gm_db": float(min(increasing, default=np.inf)),
            "gain_reduction_margin_db": float(min(decreasing)) if len(decreasing) else None,
            "pm_deg": float(min(pm)),
            "gain_boundaries_db": boundaries.tolist(),
            "phase_crossovers_rad_s": wpc.tolist(),
            "gain_crossovers_rad_s": wgc.tolist(),
        }
    return report


def synthesize(plant):
    """Set PI zeros, then normalize each return ratio at its chosen crossover."""
    data = tomllib.loads((ROOT / "data/f16/f16_sl.toml").read_text())
    elevator = data["actuators"]["elevator"]
    limits = (
        tuple(np.deg2rad(elevator["position_limits"])),
        tuple(np.deg2rad([-30, 30])), tuple(np.deg2rad([-15, 15])), (0, 1),
    )
    gains = {"rate_hz": RATE_HZ}
    for name, wc, (low, high) in zip(LOOPS, CROSSOVERS, limits):
        gains[name] = dict(kp=0.0, ki=0.0, kd=0.0, tf=1 / RATE_HZ, tt=1 / wc,
                           umin=float(low), umax=float(high))
    for i, (name, wc, zero) in enumerate(zip(LOOPS, CROSSOVERS, ZEROS)):
        gains[name].update(kp=1.0, ki=zero)
        magnitude = abs(ct.evalfr(return_ratio(plant, gains, i), 1j * wc))
        # Positive elevator is nose-down: the pitch-rate controller is negative.
        kp = (-1 if i == 0 else 1) / magnitude
        gains[name].update(kp=float(kp), ki=float(kp * zero))
    return gains


def gains_text(gains):
    lines = ["# Generated by python/control_design/design.py; SI units, angles in radians.",
             "# Outputs are absolute commands; initialize each PID with reset(trim_output).",
             f'rate_hz = {gains["rate_hz"]!r}', ""]
    for name in LOOPS:
        lines.append(f"[{name}]")
        # Ten significant digits retain design accuracy across compiler roundoff.
        for key, value in gains[name].items():
            rounded = float(format(value, ".10g"))
            lines.append(f"{key} = {rounded!r}")
        lines.append("")
    return "\n".join(lines)


def step_predictions(plant, gains):
    system = closed_loop(plant, gains)
    time = np.linspace(0, 300, 3001)  # 0.1 s overlay grid; resolves the altitude/airspeed responses.
    result = {"units": "SI perturbations from trim", "step_time_s": 0.0}
    for name, column, step in (("C1", 2, 30.48), ("C3", 3, 15.24)):
        commands = np.zeros((4, len(time)))
        commands[column] = step
        output = ct.forced_response(system, time, commands).outputs
        result[name] = dict(step_si=step, time_s=time.tolist(),
                            altitude_m=output[2].tolist(), airspeed_mps=output[3].tolist())
    return result


def generate(destination):
    destination = Path(destination)
    plant = load_plant()
    gains = tomllib.loads(gains_text(synthesize(plant)))
    report = evaluate(plant, gains)
    if not closed_loop_stable(plant, gains):
        raise ValueError("unstable closed-loop design")
    for name, loop in report.items():
        if min(loop["gm_db"], loop["gain_reduction_margin_db"] or np.inf) < 6 or loop["pm_deg"] < 45:
            raise ValueError(f"{name} fails C5: {loop}")
    if (report["pitch_rate"]["crossover_rad_s"] < 3 * report["pitch"]["crossover_rad_s"]
            or report["pitch"]["crossover_rad_s"] < 3 * report["altitude"]["crossover_rad_s"]):
        raise ValueError("nested loop crossovers are separated by less than 3x")
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "autopilot_gains.toml").write_text(gains_text(gains))
    for name, content in (("margins.json", report), ("linear_steps.json", step_predictions(plant, gains))):
        # All selected margins are finite; JSON remains portable (no Infinity/NaN).
        (destination / name).write_text(json.dumps(content, indent=2, allow_nan=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "build/control_design")
    args = parser.parse_args()
    report = generate(args.out)
    GAINS.write_bytes((args.out / "autopilot_gains.toml").read_bytes())
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
