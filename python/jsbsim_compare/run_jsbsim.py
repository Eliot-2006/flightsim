"""Run the generated JSBSim F-16 (PRD 8.6): own trim, open-loop maneuvers from trim, or init-at-state.

uv run python -m jsbsim_compare.run_jsbsim --case elev_doublet --dt 1/2400 --out out.csv
Cases: trim (hold 10 s at trim), elev_doublet, ail_pulse, rud_doublet, init (needs --state JSON).
init JSON (SI, our x13 layout [pN pE pD | vN vE vD | q0..q3 scalar-first Hamilton q_NB | p q r]):
  {"x13": [...13 SI values...], "controls": {"throttle": t, "de_deg": d, "da_deg": d, "dr_deg": d}}
Positions other than altitude (h = -pD) are ignored. Accelerations are printed (JSON, body axes, ft/s^2 and rad/s^2).
"""
import argparse
import csv
import json
import math
from fractions import Fraction
from pathlib import Path

import jsbsim

from jsbsim_compare import gen_xml

ROOT = gen_xml.REPO / "reference/jsbsim"
FT = 1 / 0.3048
# Integrators (jsbsim 1.3.3 enum: 2 Trapezoidal, 4 AdamsBashforth3). Plan 6.2: rates AB3, position rot trapezoidal, trans AB3.
INTEGRATORS = {"simulation/integrator/rate/rotational": 4, "simulation/integrator/rate/translational": 4,
               "simulation/integrator/position/rotational": 2, "simulation/integrator/position/translational": 4}
COLS = ["t_s", "V_fps", "alpha_deg", "beta_deg", "phi_deg", "theta_deg", "psi_deg", "p_dps", "q_dps", "r_dps",
        "h_ft", "pN_ft", "pE_ft", "throttle", "de_deg", "da_deg", "dr_deg"]
V0 = 502.0
T_END = 20.0
LOG_DT = 0.01  # CSV sample interval; divides both 1/2400 and 1/4800 so runs compare at identical times
# (direct-drive property, [(t_start, t_end, amplitude_deg)]) -- PRD 8.6 J3 / task spec
MANEUVERS = {
    "elev_doublet": ("aero/de-direct-deg", [(1, 2, 1.0), (2, 3, -1.0)]),
    "ail_pulse": ("aero/da-direct-deg", [(1, 2, 2.0)]),
    "rud_doublet": ("aero/dr-direct-deg", [(1, 2, 2.0), (2, 3, -2.0)]),
}


def make(dt):
    fdm = jsbsim.FGFDMExec(str(ROOT))
    fdm.set_debug_level(0)
    assert fdm.load_planet(str(ROOT / "planet/flat_nonrot.xml"), False)  # planet BEFORE model (spike)
    assert fdm.load_model("f16_sl")
    fdm.set_dt(dt)
    fdm["simulation/gravity-model"] = 0
    return fdm


def set_integrators(fdm):
    for k, v in INTEGRATORS.items():
        fdm[k] = v


def row(fdm, n0, e0, r0):
    lat, lon = fdm["position/lat-gc-rad"], fdm["position/long-gc-rad"]
    d = math.degrees
    return [fdm.get_sim_time(), fdm["velocities/vt-fps"], fdm["aero/alpha-deg"], fdm["aero/beta-deg"],
            fdm["attitude/phi-deg"], fdm["attitude/theta-deg"], fdm["attitude/psi-deg"],
            d(fdm["velocities/p-rad_sec"]), d(fdm["velocities/q-rad_sec"]), d(fdm["velocities/r-rad_sec"]),
            fdm["position/h-sl-ft"], (lat - n0) * r0, (lon - e0) * math.cos(lat) * r0, fdm["fcs/throttle-cmd-norm[0]"],
            fdm["aero/de-total-deg"], fdm["aero/da-total-deg"], fdm["aero/dr-total-deg"]]


def trim_at(fdm, h):
    fdm["ic/vt-fps"], fdm["ic/h-sl-ft"] = V0, h
    fdm["ic/lat-gc-deg"] = fdm["ic/long-gc-deg"] = fdm["ic/psi-true-deg"] = fdm["ic/gamma-deg"] = 0.0
    fdm["ic/alpha-deg"], fdm["ic/beta-deg"] = 2.5, 0.0
    fdm["propulsion/set-running"] = -1
    fdm["fcs/throttle-cmd-norm[0]"] = 0.2
    assert fdm.run_ic()
    fdm["simulation/do_simple_trim"] = 1  # full trim; raises on failure
    set_integrators(fdm)  # re-assert in case trim touched them


def euler_from_q(q0, q1, q2, q3):
    phi = math.atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2))
    theta = math.asin(max(-1.0, min(1.0, 2 * (q0 * q2 - q3 * q1))))
    psi = math.atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3))
    return map(math.degrees, (phi, theta, psi))


def init_state(fdm, spec):
    x = spec["x13"]
    c = spec["controls"]
    phi, theta, psi = euler_from_q(*x[6:10])
    fdm["ic/h-sl-ft"] = -x[2] * FT
    fdm["ic/lat-gc-deg"] = fdm["ic/long-gc-deg"] = 0.0
    fdm["ic/vn-fps"], fdm["ic/ve-fps"], fdm["ic/vd-fps"] = (v * FT for v in x[3:6])
    fdm["ic/phi-deg"], fdm["ic/theta-deg"], fdm["ic/psi-true-deg"] = phi, theta, psi
    fdm["ic/p-rad_sec"], fdm["ic/q-rad_sec"], fdm["ic/r-rad_sec"] = x[10:13]
    fdm["propulsion/set-running"] = -1
    fdm["fcs/throttle-cmd-norm[0]"] = c["throttle"]
    fdm["aero/de-direct-deg"], fdm["aero/da-direct-deg"], fdm["aero/dr-direct-deg"] = c["de_deg"], c["da_deg"], c["dr_deg"]
    assert fdm.run_ic()
    fdm["fcs/throttle-cmd-norm[0]"] = c["throttle"]
    keys = {"udot_fps2": "accelerations/udot-ft_sec2", "vdot_fps2": "accelerations/vdot-ft_sec2", "wdot_fps2": "accelerations/wdot-ft_sec2",
            "pdot_rps2": "accelerations/pdot-rad_sec2", "qdot_rps2": "accelerations/qdot-rad_sec2", "rdot_rps2": "accelerations/rdot-rad_sec2"}
    return {k: fdm[v] for k, v in keys.items()}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=["trim", *MANEUVERS, "init"])
    ap.add_argument("--dt", default="1/2400", help="step, e.g. 1/2400")
    ap.add_argument("--h-ft", type=float, default=0.0, help="trim altitude (ft); spec uses 0 and 10000")
    ap.add_argument("--state", help="JSON for --case init")
    ap.add_argument("--out", help="CSV path (metadata JSON written beside it)")
    a = ap.parse_args(argv)
    dt = float(Fraction(a.dt))
    fdm = make(dt)
    set_integrators(fdm)
    meta = {"jsbsim_version": jsbsim.__version__, "case": a.case, "dt": a.dt, "log_dt_s": LOG_DT, "integrators": INTEGRATORS,
            "gravity_model": 0, "planet": "reference/jsbsim/planet/flat_nonrot.xml", "xcg_cbar": 0.35,
            "pN_pE": f"arc length from start at radius R+h0, R={gen_xml.R_FT} ft"}
    if a.case == "init":
        acc = init_state(fdm, json.loads(Path(a.state).read_text()))
        print(json.dumps(acc, indent=1))
        meta["accelerations"] = acc
        if a.out:
            Path(a.out).with_suffix(".json").write_text(json.dumps(meta, indent=1))
        return acc
    trim_at(fdm, a.h_ft)
    meta["trim"] = {k: fdm[v] for k, v in {"alpha_deg": "aero/alpha-deg", "theta_deg": "attitude/theta-deg", "throttle": "fcs/throttle-cmd-norm[0]",
                                           "de_deg": "aero/de-total-deg", "da_deg": "aero/da-total-deg", "dr_deg": "aero/dr-total-deg",
                                           "udot": "accelerations/udot-ft_sec2", "wdot": "accelerations/wdot-ft_sec2", "qdot": "accelerations/qdot-rad_sec2"}.items()}
    meta["trim"]["h_ft"] = a.h_ft
    meta["integrator_readback"] = {k: fdm[k] for k in INTEGRATORS}
    r0 = gen_xml.R_FT + fdm["position/h-sl-ft"]
    n0, e0 = fdm["position/lat-gc-rad"], fdm["position/long-gc-rad"]
    t0 = fdm.get_sim_time()
    steps, every = round(T_END / dt), round(LOG_DT / dt)
    prop, windows = MANEUVERS.get(a.case, (None, []))
    rows = [row(fdm, n0, e0, r0)]
    rows[0][0] -= t0
    for i in range(1, steps + 1):
        t = i * dt
        if prop:
            fdm[prop] = sum(amp for s, e, amp in windows if s <= t < e)
        assert fdm.run()
        if i % every == 0:
            r = row(fdm, n0, e0, r0)
            r[0] = t
            rows.append(r)
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(COLS)
            w.writerows(rows)
        Path(a.out).with_suffix(".json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta["trim"], indent=1))


if __name__ == "__main__":
    main()
