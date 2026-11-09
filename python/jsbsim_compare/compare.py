"""JSBSim matched cross-validation: J1 point check, J2 trim, J3 open-loop responses, discrepancy budget (PRD 8.6).

    uv run python -m jsbsim_compare.compare [--out out/jsbsim] [--bin build/release/flightsim] [--check] [--regen]

Our side is the C++ CLI (`flightsim run|trim` on scenarios/jsbsim/*.toml); the JSBSim side is the committed
reference/jsbsim/results/*.csv (J3, trim) plus a live JSBSim init-at-state run (J1). Python never integrates our model.
--check regenerates the JSBSim XML and every JSBSim CSV and asserts they equal the committed ones, then asserts the
J1-J3 criteria. --regen rewrites the committed CSVs. Writes <out>/compare.json (consumed by the report).
"""

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

import numpy as np

from jsbsim_compare import gen_xml, run_jsbsim

REPO = gen_xml.REPO
RESULTS = REPO / "reference/jsbsim/results"
SCEN = REPO / "scenarios/jsbsim"
FT = 1 / 0.3048
H_J3_FT = 10000  # J3 altitude: from h = 0 the maneuvers leave the thrust table's lower altitude edge (see scenarios)

# ---- criteria. Origins: PRD 8.6 (J1-J3) -------------------------------------------------------------------------
J1_UDOT_WDOT_FPS2 = 1e-3  # PRD 8.6 J1 (initial target; exceeding it blocks until root-caused)
J1_QDOT_RPS2 = 1e-5  # PRD 8.6 J1
J2_REL = 1e-3  # PRD 8.6 J2: alpha, de, throttle within 0.1 %
J3_REL = 0.02  # PRD 8.6 J3: max difference <= 2 % of each signal's peak deviation
NOT_EXCITED = 1e-6  # a signal whose JSBSim peak deviation is below this (deg, deg/s, ft, ft/s) is round-off, not a response
CSV_REL = 1e-9  # plan Task 6.3: regenerated JSBSim CSVs match the committed ones within 1e-9 relative (to column peak)

CASES = {  # name -> (label, signals shown in the overlay figure)
    "elev_doublet": ("Elevator doublet ±1°, 1–2 s / 2–3 s", ["V", "alpha", "q", "theta"]),
    "ail_pulse": ("Aileron pulse 2°, 1–2 s", ["beta", "p", "r", "phi"]),
    "rud_doublet": ("Rudder doublet ±2°, 1–2 s / 2–3 s", ["beta", "p", "r", "phi"]),
}
# signal -> (CSV column, label, unit)
SIGNALS = {
    "V": ("V_fps", "ΔV", "ft/s"), "alpha": ("alpha_deg", "Δα", "deg"), "q": ("q_dps", "Δq", "deg/s"),
    "theta": ("theta_deg", "Δθ", "deg"), "h": ("h_ft", "Δh", "ft"), "beta": ("beta_deg", "Δβ", "deg"),
    "p": ("p_dps", "Δp", "deg/s"), "r": ("r_dps", "Δr", "deg/s"), "phi": ("phi_deg", "Δφ", "deg"),
}


# ---------------------------------------------------------------------------------------------- helpers
def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def flightsim(binary, *args):
    r = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{binary} {' '.join(map(str, args))} failed:\n{r.stdout}{r.stderr}")
    return r.stdout


def our_run(binary, scenario, out):
    """Run a scenario through `flightsim run`; returns (csv columns, metadata)."""
    out = Path(out)
    flightsim(binary, "run", scenario, "--out", out)
    meta = json.loads((out / "metadata.json").read_text())
    assert not meta["out_of_domain"], f"{scenario}: our run left a table domain: {meta['first_violation']}"
    return read_csv(out / "trajectory.csv"), meta


def our_trim(binary, scenario, out):
    flightsim(binary, "trim", scenario, "--out", out)
    d = json.loads(Path(out).read_text())
    assert d["trim"]["ok"] and not d["trim"]["lookup_out_of_domain"], f"{scenario}: trim failed"
    return d["trim"], d


def deviation(cols, key):
    return cols[key] - cols[key][0]


def metric(ours, jsb, key, hold=None):
    """max|d_ours - d_jsb| / max|d_jsb| for one signal; None if JSBSim shows no response. `hold` = JSBSim no-input drift
    to subtract from d_jsb (curvature correction)."""
    do, dj = deviation(ours, key), deviation(jsb, key)
    if hold is not None:
        dj = dj - deviation(hold, key)
    peak = np.abs(dj).max()
    return None if peak < NOT_EXCITED else float(np.abs(do - dj).max() / peak)


def metrics_for(ours, jsb, hold=None):
    return {s: metric(ours, jsb, SIGNALS[s][0], hold) for s in SIGNALS}


def worst(m):
    vals = [v for v in m.values() if v is not None]
    return max(vals) if vals else 0.0


# ---------------------------------------------------------------------------------------------- J1 / J2
def j1(binary, work):
    out = []
    for h in (0, 10000):
        tr, _ = our_trim(binary, SCEN / f"trim_h{h}.toml", work / f"trim_h{h}.json")
        state = work / f"state_h{h}.json"
        state.write_text(json.dumps({"x13": tr["x17"][:13], "controls": {
            "throttle": tr["throttle"], "de_deg": tr["elevator_deg"], "da_deg": 0.0, "dr_deg": 0.0}}))
        acc = run_jsbsim.main(["--case", "init", "--state", str(state)])
        ours = tr["accel_body"]
        # Decomposition: JSBSim's density differs from ours by (rho_j - rho_o)/rho_o, so its lift (~weight) differs by that
        # fraction and w-dot (down positive) by -g * fraction. Reported next to the measured value.
        jatm = _jsbsim_atmosphere(h)
        rel_rho = (jatm["rho"] - tr["atmosphere"]["rho"]) / tr["atmosphere"]["rho"]
        g_fps2 = 9.80665 * (gen_xml.R_FT / FT / (gen_xml.R_FT / FT + h / FT)) ** 2 * FT
        out.append({
            "rho_rel_diff": rel_rho, "wdot_from_rho_fps2": -g_fps2 * rel_rho,
            "h_ft": h,
            "jsbsim": {"udot_fps2": acc["udot_fps2"], "wdot_fps2": acc["wdot_fps2"], "qdot_rps2": acc["qdot_rps2"],
                       "vdot_fps2": acc["vdot_fps2"], "pdot_rps2": acc["pdot_rps2"], "rdot_rps2": acc["rdot_rps2"]},
            "ours": {"udot_fps2": ours["udot_mps2"] * FT, "wdot_fps2": ours["wdot_mps2"] * FT, "qdot_rps2": ours["qdot_rps2"]},
            "pass": abs(acc["udot_fps2"]) < J1_UDOT_WDOT_FPS2 and abs(acc["wdot_fps2"]) < J1_UDOT_WDOT_FPS2
            and abs(acc["qdot_rps2"]) < J1_QDOT_RPS2,
        })
    return out


def j2(binary, work):
    rows = []
    for h in (0, 10000):
        tr, _ = our_trim(binary, SCEN / f"trim_h{h}.toml", work / f"trim_h{h}.json")
        jt = json.loads((RESULTS / f"trim_h{h}_dt2400.json").read_text())["trim"]
        pairs = {"alpha_deg": (tr["alpha_deg"], jt["alpha_deg"]), "de_deg": (tr["elevator_deg"], jt["de_deg"]),
                 "throttle": (tr["throttle"], jt["throttle"])}
        r = {"h_ft": h}
        for k, (a, b) in pairs.items():
            r[k] = {"ours": a, "jsbsim": b, "rel": abs(a - b) / abs(b)}
        r["pass"] = all(r[k]["rel"] <= J2_REL for k in pairs)
        rows.append(r)
    return rows


# ---------------------------------------------------------------------------------------------- budget variants
def variant_text(base: str, kind: str, jsb_trim=None):
    """Our scenario text with ONE factor changed. `base` is a scenarios/jsbsim/*.toml text."""
    t = base.replace('data = "data/f16/f16_sl.toml"', f'data = "{REPO / "data/f16/f16_sl.toml"}"')
    if kind == "gravity_constant":  # drop the inverse-square law (g0 held)
        t = re.sub(r"R = [0-9.e+]+ \}", "R = inf }", t)
    elif kind == "atmosphere_sl_adc":  # sensitivity probe: a different atmosphere, scaled to the measured rho difference
        t = t.replace('atmosphere = "us76"', 'atmosphere = "sl_adc"')
    elif kind == "jsbsim_trim_state":  # start from JSBSim's own trim (alpha, theta, throttle, de) instead of ours
        g = tomllib.loads((REPO / "data/f16/f16_sl.toml").read_text())["engine"]["throttle_gearing"]
        thr = jsb_trim["throttle"]
        power = g["low_slope"] * thr if thr <= g["breakpoint"] else g["high_slope"] * thr + g["high_intercept"]
        x13 = [502 / FT, math.radians(jsb_trim["alpha_deg"]), 0, 0, math.radians(jsb_trim["theta_deg"]), 0, 0, 0, 0, 0, 0,
               H_J3_FT / FT, power]
        t = t.replace("[initial]\nfrom_trim = true", "[initial]\nfrom_trim = false\nx13 = [" + ", ".join(repr(x) for x in x13)
                      + f"]\nthrottle = {thr!r}\nde_deg = {jsb_trim['de_deg']!r}\nda_deg = 0.0\ndr_deg = 0.0")
        assert "from_trim = false" in t
    else:
        raise ValueError(kind)
    return t


def analytic_curvature(V_mps, R_m, h_m, t):
    """Kinematic spherical-Earth effect, valid while aerodynamic feedback is still small (first seconds): JSBSim's body holds
    its inertial attitude (q = 0, planet not rotating) and flies straight, so relative to the curving local horizontal it
    climbs: dh = V^2 t^2 / (2 (R+h)), dtheta = +V t / (R+h), and gravity (pointing at the centre) decelerates it along track:
    dV = -g V t^2 / (2 (R+h)). Our flat-Earth model has none of this; compared against JSBSim's measured no-input drift."""
    r = R_m + h_m
    g = 9.80665 * (R_m / r) ** 2
    return {"h_ft": (V_mps * t) ** 2 / (2 * r) * FT, "theta_deg": math.degrees(V_mps * t / r), "V_fps": -g * V_mps * t**2 / (2 * r) * FT}


# ---------------------------------------------------------------------------------------------- main computation
def compute(binary, out):
    out = Path(out)
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    res = {"criteria": {"J1_udot_wdot_fps2": J1_UDOT_WDOT_FPS2, "J1_qdot_rps2": J1_QDOT_RPS2, "J2_rel": J2_REL, "J3_rel": J3_REL,
                        "not_excited_below": NOT_EXCITED}}
    res["J1"], res["J2"] = j1(binary, work), j2(binary, work)
    jt = json.loads((RESULTS / f"trim_h{H_J3_FT}_dt2400.json").read_text())["trim"]
    hold = {n: read_csv(RESULTS / f"trim_h{H_J3_FT}_dt{n}.csv") for n in (2400, 4800)}
    trim10, trim10_full = our_trim(binary, SCEN / f"trim_h{H_J3_FT}.toml", work / "trim_h10000.json")

    j3, series, budget = {}, {}, {k: {} for k in (
        "jsbsim_dt_halving", "our_dt_halving", "gravity_constant", "atmosphere", "curvature", "residual_trim_mismatch",
        "sampling_alignment", "table_edges", "all_explained")}
    atm = {}
    for case in CASES:
        stem = f"{case}_h{H_J3_FT}"
        base_text = (SCEN / f"{stem}_dt2400.toml").read_text()
        ours, meta = our_run(binary, SCEN / f"{stem}_dt2400.toml", work / f"{case}_ours")
        ours_h, _ = our_run(binary, SCEN / f"{stem}_dt1200.toml", work / f"{case}_ours_dt1200")
        jsb, jsb_h = read_csv(RESULTS / f"{stem}_dt2400.csv"), read_csv(RESULTS / f"{stem}_dt4800.csv")
        dt_err = np.abs(ours["t_s"] - jsb["t_s"]).max()
        assert dt_err < 1e-9, f"{case}: sample times differ by {dt_err}"
        assert np.abs(jsb["t_s"] - jsb_h["t_s"]).max() < 1e-9
        base = metrics_for(ours, jsb)
        j3[case] = {"metrics": base, "worst": worst(base), "pass": worst(base) <= J3_REL,
                    "peak": {s: float(np.abs(deviation(jsb, SIGNALS[s][0])).max()) for s in SIGNALS}}
        series[case] = {"t_s": ours["t_s"].tolist(), **{s: {"ours": deviation(ours, SIGNALS[s][0]).tolist(),
                                                           "jsbsim": deviation(jsb, SIGNALS[s][0]).tolist()} for s in SIGNALS}}

        def effect(m):  # change of each J3 metric when one factor is toggled
            return {s: (None if m[s] is None or base[s] is None else m[s] - base[s]) for s in SIGNALS}

        # JSBSim discretization: its own dt/2 run, compared both directly and through the J3 metric
        budget["jsbsim_dt_halving"][case] = {"effect": effect(metrics_for(ours, jsb_h)),
                                             "self": {s: _self(jsb, jsb_h, SIGNALS[s][0]) for s in SIGNALS}}
        budget["our_dt_halving"][case] = {"effect": effect(metrics_for(ours_h, jsb)),
                                          "self": {s: _self(ours, ours_h, SIGNALS[s][0]) for s in SIGNALS}}
        # gravity: constant g in ours
        v = work / f"{case}_gravity_constant.toml"
        v.write_text(variant_text(base_text, "gravity_constant"))
        g_cols, _ = our_run(binary, v, work / f"{case}_gravity_constant")
        budget["gravity_constant"][case] = {"effect": effect(metrics_for(g_cols, jsb))}
        # atmosphere: sensitivity probe with sl_adc, scaled to the measured US76 rho mismatch (below)
        v = work / f"{case}_atmosphere_sl_adc.toml"
        v.write_text(variant_text(base_text, "atmosphere_sl_adc"))
        a_cols, _ = our_run(binary, v, work / f"{case}_atmosphere_sl_adc")
        budget["atmosphere"][case] = {"probe_effect": effect(metrics_for(a_cols, jsb))}
        # curvature: JSBSim's own no-input drift (spherical Earth + its trim residual) subtracted from its maneuver
        budget["curvature"][case] = {"effect": effect(metrics_for(ours, jsb, hold[2400]))}
        # residual trim mismatch: start ours from JSBSim's trim values
        v = work / f"{case}_jsbsim_trim_state.toml"
        v.write_text(variant_text(base_text, "jsbsim_trim_state", jt))
        t_cols, _ = our_run(binary, v, work / f"{case}_jsbsim_trim_state")
        budget["residual_trim_mismatch"][case] = {"effect": effect(metrics_for(t_cols, jsb))}
        # everything above explained at once: JSBSim trim state, its no-input drift removed, JSBSim at the smaller step
        budget["all_explained"][case] = {"metrics": metrics_for(t_cols, jsb_h, hold[4800])}
        # sampling: both sides log on the same 0.01 s grid; bound = max|d(delta)/dt| * time mismatch / peak
        rate = {s: float(np.abs(np.gradient(deviation(jsb, SIGNALS[s][0]), jsb["t_s"])).max()) for s in SIGNALS}
        budget["sampling_alignment"][case] = {"time_mismatch_s": float(dt_err), "effect": {
            s: (None if j3[case]["peak"][s] < NOT_EXCITED else rate[s] * float(dt_err) / j3[case]["peak"][s]) for s in SIGNALS}}
        # table edges: our monitor is asserted false in our_run; check JSBSim stayed inside the same axes
        # (effect is exactly zero: with both sides inside the table domain, clamp-vs-extrapolate never acts)
        budget["table_edges"][case] = {"domain": _jsbsim_domain(jsb), "effect": {
            s: (None if j3[case]["peak"][s] < NOT_EXCITED else 0.0) for s in SIGNALS}}

    # atmosphere mismatch at the J3 altitude, from JSBSim's atmosphere/* properties vs our us76 (via trim.json)
    jatm = _jsbsim_atmosphere(H_J3_FT)
    our_atm = trim10_full["trim"]["atmosphere"]
    sl_atm = _our_sladc(binary, work)
    atm = {"h_ft": H_J3_FT, "ours": our_atm, "jsbsim": jatm,
           "rel": {k: (our_atm[k] - jatm[k]) / jatm[k] for k in ("rho", "T", "a")},
           "probe_rel_rho": (sl_atm["rho"] - our_atm["rho"]) / our_atm["rho"]}
    for case in CASES:  # scale the sl_adc probe to the actual mismatch (metrics respond linearly to a small rho change)
        s_ = abs(atm["rel"]["rho"] / atm["probe_rel_rho"])
        budget["atmosphere"][case]["effect"] = {k: (None if v is None else v * s_)
                                                for k, v in budget["atmosphere"][case]["probe_effect"].items()}
        budget["atmosphere"][case]["scale"] = s_
    res["atmosphere"] = atm
    res["J3"], res["budget"], res["series"] = j3, budget, series
    # curvature: kinematic estimate vs the measured JSBSim no-input drift
    V = 502 / FT
    h = hold[2400]
    res["curvature"] = {"t_s": [5.0, 20.0], "analytic": [analytic_curvature(V, gen_xml.R_FT / FT, H_J3_FT / FT, t) for t in (5.0, 20.0)],
                        "jsbsim_hold": [{k: float(h[k][i] - h[k][0]) for k in ("h_ft", "theta_deg", "V_fps")} for i in (500, 2000)]}
    res["jsbsim_trim"] = jt
    res["our_trim"] = trim10
    return res


def _self(a, b, key):
    """max|d_a - d_b| / max|d_a|: how much one sim's own deviation moves when its step is halved."""
    da, db = deviation(a, key), deviation(b, key)
    peak = np.abs(da).max()
    return None if peak < NOT_EXCITED else float(np.abs(da - db).max() / peak)


def _jsbsim_domain(jsb):
    d = tomllib.loads((REPO / "data/f16/f16_sl.toml").read_text())["aero"]
    a_ax, b_ax, e_ax = d["cx"]["alpha"], d["dlda"]["beta"], d["cx"]["elevator"]
    rng = {"alpha_deg": (min(a_ax), max(a_ax), jsb["alpha_deg"]), "abs_beta_deg": (0, max(abs(x) for x in b_ax), np.abs(jsb["beta_deg"])),
           "de_deg": (min(e_ax), max(e_ax), jsb["de_deg"])}
    out = {k: {"min": float(v.min()), "max": float(v.max()), "table": [lo, hi], "inside": bool(v.min() >= lo and v.max() <= hi)}
           for k, (lo, hi, v) in rng.items()}
    assert all(r["inside"] for r in out.values()), f"JSBSim run left a table domain: {out}"
    return out


def _jsbsim_atmosphere(h_ft):
    fdm = run_jsbsim.make(1 / 2400)
    fdm["ic/vt-fps"], fdm["ic/h-sl-ft"] = run_jsbsim.V0, h_ft
    fdm["ic/lat-gc-deg"] = fdm["ic/long-gc-deg"] = fdm["ic/psi-true-deg"] = fdm["ic/gamma-deg"] = 0.0
    assert fdm.run_ic()
    # JSBSim: slug/ft^3, degR, ft/s -> SI
    return {"rho": fdm["atmosphere/rho-slugs_ft3"] * 14.59390294 / 0.3048**3, "T": fdm["atmosphere/T-R"] / 1.8,
            "a": fdm["atmosphere/a-fps"] * 0.3048}


def _our_sladc(binary, work):
    t = variant_text((SCEN / f"trim_h{H_J3_FT}.toml").read_text(), "atmosphere_sl_adc")
    (work / "trim_sladc.toml").write_text(t)
    return our_trim(binary, work / "trim_sladc.toml", work / "trim_sladc.json")[1]["trim"]["atmosphere"]


# ---------------------------------------------------------------------------------------------- check / regen
def reference_runs():
    """(case, dt) -> run_jsbsim argv for every committed JSBSim result."""
    for dt in (2400, 4800):
        for h in (0, 10000):
            yield f"trim_h{h}_dt{dt}", ["--case", "trim", "--dt", f"1/{dt}", "--h-ft", str(h)]
        for case in CASES:
            yield f"{case}_h{H_J3_FT}_dt{dt}", ["--case", case, "--dt", f"1/{dt}", "--h-ft", str(H_J3_FT)]


def regenerate_jsbsim(dest):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    for name, argv in reference_runs():
        run_jsbsim.main([*argv, "--out", str(dest / f"{name}.csv")])


# Physical floor for each column's scale. A column whose peak is below its floor carries no response: in the
# symmetric cases the lateral channels sit at round-off (~1e-15), which differs between platforms' libm.
# Relative-to-peak would turn that noise into O(1) "differences". Floors are far below any excited J3 peak.
CSV_FLOOR = {"ft": 1e-3, "fps": 1e-6, "deg": 1e-6, "dps": 1e-6}


def _floor(col):
    return next((v for suffix, v in CSV_FLOOR.items() if col.endswith("_" + suffix)), 1e-9)


def csv_matches(a, b):
    """Regenerated vs committed CSV: same columns; every element within CSV_REL of max(column peak, physical floor).

    Angles are compared modulo 360 deg (psi wraps)."""
    ca, cb = read_csv(a), read_csv(b)
    assert ca.keys() == cb.keys(), f"{a.name}: columns differ"
    bad = []
    for k in ca:
        diff = ca[k] - cb[k]
        if k.endswith("_deg"):
            diff = (diff + 180.0) % 360.0 - 180.0
        scale = max(np.abs(cb[k]).max(), _floor(k))
        err = float(np.abs(diff).max() / scale)
        if not err <= CSV_REL:
            bad.append((k, err))
    return bad


def check(res, tmp):
    fails = []
    for r in res["J1"]:
        if not r["pass"]:
            fails.append(f"J1 h={r['h_ft']} ft: {r['jsbsim']}")
    for r in res["J2"]:
        if not r["pass"]:
            fails.append(f"J2 h={r['h_ft']} ft: {r}")
    for case, r in res["J3"].items():
        if not r["pass"]:
            fails.append(f"J3 {case}: worst {r['worst']:.4g} > {J3_REL}: {r['metrics']}")
    # generated XML is exactly the committed XML
    gen_xml.generate(tmp / "xml")
    for rel in ("aircraft/f16_sl/f16_sl.xml", "engine/dummy.xml", "engine/direct.xml", "planet/flat_nonrot.xml"):
        if (tmp / "xml" / rel).read_text() != (gen_xml.REPO / "reference/jsbsim" / rel).read_text():
            fails.append(f"regenerated {rel} differs from the committed file")
    regenerate_jsbsim(tmp / "csv")
    for name, _ in reference_runs():
        bad = csv_matches(tmp / "csv" / f"{name}.csv", RESULTS / f"{name}.csv")
        if bad:
            fails.append(f"{name}.csv differs from the committed CSV beyond {CSV_REL}: {bad}")
    return fails


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(REPO / "out/jsbsim"))
    ap.add_argument("--bin", default=os.environ.get("FLIGHTSIM_BIN", str(REPO / "build/release/flightsim")))
    ap.add_argument("--check", action="store_true", help="also regenerate JSBSim XML+CSVs and assert J1-J3 and CSV equality")
    ap.add_argument("--regen", action="store_true", help="rewrite the committed reference/jsbsim CSVs from JSBSim")
    a = ap.parse_args(argv)
    if a.regen:
        gen_xml.generate(gen_xml.REPO / "reference/jsbsim")
        regenerate_jsbsim(RESULTS)
    res = compute(a.bin, a.out)
    Path(a.out).mkdir(parents=True, exist_ok=True)
    (Path(a.out) / "compare.json").write_text(json.dumps(res))
    print_summary(res)
    if a.check:
        with tempfile.TemporaryDirectory() as t:
            fails = check(res, Path(t))
        if fails:
            print("\nCHECK FAILED:\n  " + "\n  ".join(fails))
            return 1
        print("\nCHECK PASSED: XML regenerates identically, JSBSim CSVs match, J1-J3 criteria met")
    return 0


def print_summary(res):
    for r in res["J1"]:
        j = r["jsbsim"]
        print(f"J1 h={r['h_ft']:>5} ft  JSBSim udot={j['udot_fps2']:.3e} wdot={j['wdot_fps2']:.3e} ft/s2 qdot={j['qdot_rps2']:.3e} rad/s2  "
              f"{'PASS' if r['pass'] else 'FAIL'}")
    for r in res["J2"]:
        print(f"J2 h={r['h_ft']:>5} ft  " + "  ".join(f"{k} rel={r[k]['rel']:.2e}" for k in ("alpha_deg", "de_deg", "throttle"))
              + f"  {'PASS' if r['pass'] else 'FAIL'}")
    for case, r in res["J3"].items():
        print(f"J3 {case:13s} worst={r['worst']:.3e}  " + " ".join(f"{s}={'-' if v is None else format(v, '.1e')}" for s, v in r["metrics"].items())
              + f"  {'PASS' if r['pass'] else 'FAIL'}")
    for src, cases in res["budget"].items():
        w = 0.0
        for c in cases.values():
            e = c.get("effect") or c.get("metrics") or {}
            w = max([w] + [abs(v) for v in e.values() if v is not None])
        print(f"budget {src:24s} worst |effect on J3 metric| = {w:.2e}")


if __name__ == "__main__":
    sys.exit(main())
