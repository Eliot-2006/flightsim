"""J4 (PRD 8.6 appendix): JSBSim's stock `f16` (from the pinned jsbsim wheel) trimmed at the matched condition.

    uv run python -m jsbsim_compare.stock_f16 [--out out/jsbsim] [--bin build/release/flightsim]

Trim only. The stock flight-control system overwrites direct surface writes (spikes/jsbsim/FINDINGS.md Q7), so no
open-loop response is compared. The stock aircraft is a different model (General Dynamics F-16A, Erik Hofman, GPL;
NASA TP-1538 + Dash-1 + Murray data, not the S&L reduced tables), so the numbers are a model-form difference with no
pass/fail. Ours is the matched model of J2 (scenarios/jsbsim/trim_h*.toml: US76, inverse-square gravity, no rotor,
S&L x_cg = 0.35). Stock keeps JSBSim's default planet (rotating WGS84, gravity model 1): J4 is not matched on Earth.
Writes <out>/stock_f16.json (consumed by the report).
"""

import argparse
import json
import math
import os
import tempfile
from pathlib import Path

import jsbsim

from jsbsim_compare import compare, gen_xml

REPO = gen_xml.REPO
SCEN = REPO / "scenarios/jsbsim"
V0_FPS = 502.0
LBF_PER_N = 1 / 4.4482216152605  # data/f16 conversion (PRD 4.3)
# Ground contact would stop FGTrim at sea level (the file has gear and structure contact points): the terrain is moved
# 5,000 ft below the aircraft. It does not change the atmosphere, which is a function of MSL altitude.
TERRAIN_FT = -5000.0
# Properties read back after trim: report key -> JSBSim property
READ = {
    "alpha_deg": "aero/alpha-deg", "theta_deg": "attitude/theta-deg", "mach": "velocities/mach",
    "V_fps": "velocities/vt-fps", "qbar_psf": "aero/qbar-psf",
    "elevator_rad": "fcs/elevator-pos-rad", "elevator_deg": "fcs/elevator-pos-deg",
    "elevator_cmd_norm": "fcs/elevator-cmd-norm", "pitch_trim_cmd_norm": "fcs/pitch-trim-cmd-norm",
    "throttle_cmd_norm": "fcs/throttle-cmd-norm[0]", "throttle_pos_norm": "fcs/throttle-pos-norm[0]",
    "thrust_lbf": "propulsion/engine[0]/thrust-lbs", "weight_lbf": "inertia/weight-lbs",
    "cg_x_in": "inertia/cg-x-in", "aerorp_x_in": "metrics/aero-rp-x-in", "cbar_ft": "metrics/cbarw-ft",
    "fuel_lbf": "propulsion/total-fuel-lbs",
    "gear_pos_norm": "gear/gear-pos-norm", "CDgear": "aero/coefficient/CDgear",
    "udot_fps2": "accelerations/udot-ft_sec2", "wdot_fps2": "accelerations/wdot-ft_sec2",
    "qdot_rps2": "accelerations/qdot-rad_sec2",
}


def stock_trim(h_ft):
    """Full JSBSim trim (do_simple_trim = 1) of the stock f16 at 502 ft/s, level, gear up, default planet."""
    fdm = jsbsim.FGFDMExec(None)  # the wheel's own aircraft/ and engine/ directories
    fdm.set_debug_level(0)
    assert fdm.load_model("f16")
    for _ in range(2):  # the gear takes 5 s of FCS time to retract: trim, then re-initialise and trim from retracted
        fdm["ic/terrain-elevation-ft"] = TERRAIN_FT
        fdm["ic/vt-fps"], fdm["ic/h-sl-ft"] = V0_FPS, h_ft
        fdm["ic/gamma-deg"] = fdm["ic/beta-deg"] = fdm["ic/psi-true-deg"] = 0.0
        fdm["ic/alpha-deg"] = 3.0
        fdm["propulsion/set-running"] = -1
        fdm["fcs/throttle-cmd-norm[0]"] = 0.3
        fdm["gear/gear-cmd-norm"] = 0
        assert fdm.run_ic()
        fdm["simulation/do_simple_trim"] = 1  # full trim; raises TrimFailureError on failure
    out = {k: fdm[p] for k, p in READ.items()}
    out["h_ft"] = h_ft
    out["gravity_model"] = fdm["simulation/gravity-model"]
    out["jsbsim_version"] = jsbsim.__version__
    return out


def our_trim_at(binary, h_ft, work):
    """Our matched-model trim plus its thrust, which only the run CSV carries (one step from the trim state)."""
    tr, _ = compare.our_trim(binary, SCEN / f"trim_h{h_ft}.toml", work / f"trim_h{h_ft}.json")
    scenario = work / f"trim_h{h_ft}_run.toml"
    scenario.write_text((SCEN / f"trim_h{h_ft}.toml").read_text().replace(
        "data = \"data/f16/f16_sl.toml\"", f"data = \"{REPO / 'data/f16/f16_sl.toml'}\"")
        + "\n[sim]\nintegrator = \"rk4\"\ndt_s = 0.01\nt_final_s = 0.01\noutput_every = 1\n")
    # Not compare.our_run: one RK4 step from a sea-level trim drops h by ~1e-19 m, which the engine table (domain
    # h >= 0) flags. Only the t = 0 row is used, and trim itself asserts no violation.
    compare.flightsim(binary, "run", scenario, "--out", work / f"trim_h{h_ft}_run")
    cols = compare.read_csv(work / f"trim_h{h_ft}_run" / "trajectory.csv")
    return {"alpha_deg": tr["alpha_deg"], "elevator_deg": tr["elevator_deg"], "throttle": tr["throttle"],
            "thrust_lbf": cols["thrust_N"][0] * LBF_PER_N}


def compute(binary, work):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    rows = []
    for h in (0, 10000):
        s, o = stock_trim(h), our_trim_at(binary, h, work)
        stock = {"alpha_deg": s["alpha_deg"], "elevator_deg": s["elevator_deg"], "throttle": s["throttle_cmd_norm"],
                 "thrust_lbf": s["thrust_lbf"]}
        rows.append({"h_ft": h, "stock": stock, "ours": o, "stock_full": s,
                     "difference": {k: stock[k] - o[k] for k in o}})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(REPO / "out/jsbsim"))
    ap.add_argument("--bin", default=os.environ.get("FLIGHTSIM_BIN", str(REPO / "build/release/flightsim")))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as t:
        rows = compute(a.bin, Path(t))
    (out / "stock_f16.json").write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(f"J4 h={r['h_ft']:>5} ft  " + "  ".join(
            f"{k}: stock {r['stock'][k]:.4g} ours {r['ours'][k]:.4g}" for k in r["ours"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
