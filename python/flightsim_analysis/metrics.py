"""Collect every number the report quotes into one flat dict (out/metrics.json). Reads CLI outputs only;
the arithmetic here is post-processing of CSV columns (ratios, maxima, counts), not physics."""

import csv
import json
import math
from pathlib import Path

INTEGRATORS = ("euler", "rk4", "symplectic_euler", "lie_euler", "rkmk4")
STUDY = "core_convergence"


def read_csv(path: Path) -> dict[str, list[float]]:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def sign_changes(xs) -> int:
    signs = [x > 0 for x in xs if x != 0]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def max_rel_drift(xs) -> float:
    return max(abs(x / xs[0] - 1) for x in xs)


def max_vec_rel_drift(vs) -> float:
    n0 = math.hypot(*vs[0])
    return max(math.hypot(*(a - b for a, b in zip(v, vs[0]))) / n0 for v in vs)


def trajectory(out: Path, name: str) -> dict[str, list[float]]:
    return read_csv(out / "runs" / name / "trajectory.csv")


def hn(tr):
    return list(zip(tr["HN_x_kgm2ps"], tr["HN_y_kgm2ps"], tr["HN_z_kgm2ps"]))


def build(out: Path) -> dict:
    m: dict = {}
    meta = json.loads((out / "runs" / "rb_dzhanibekov" / "metadata.json").read_text())
    m["env.version"], m["env.git_sha"] = meta["flightsim_version"], meta["git_sha"]
    m["env.compiler"], m["env.build_type"] = meta["compiler"], meta["build_type"]

    # --- asserted convergence orders (orders.json) ---
    orders = json.loads((out / "studies" / STUDY / "orders.json").read_text())
    m["orders.all_passed"] = orders["all_assertions_passed"]
    rows, first, fourth, n_in = [], [], [], 0
    for o in orders["orders"]:
        for q in o["metrics"]:
            if not q["asserted"]:
                continue
            key = f"order.{o['case']}.{o['integrator']}.{q['metric']}"
            m[key], m[key + ".lo"], m[key + ".hi"] = q["order"], q["lower_bound"], q["upper_bound"]
            (first if q["upper_bound"] < 2 else fourth).append(q["order"])
            m["orders.first.lo" if q["upper_bound"] < 2 else "orders.fourth.lo"] = q["lower_bound"]
            m["orders.first.hi" if q["upper_bound"] < 2 else "orders.fourth.hi"] = q["upper_bound"]
            n_in += q["within_bounds"]
            note = " (unprojected q)" if o["case"] == "R1" and not o["renormalize"] else ""
            rows.append(f"| {o['case']} | `{o['integrator']}` | {q['metric']}{note} | {q['order']:.3f} | "
                        f"[{q['lower_bound']}, {q['upper_bound']}] | {'pass' if q['within_bounds'] else 'FAIL'} |")
    m["orders.asserted"], m["orders.in_bounds"] = len(rows), n_in
    m["orders.first.min"], m["orders.first.max"] = min(first), max(first)
    m["orders.fourth.min"], m["orders.fourth.max"] = min(fourth), max(fourth)
    m["orders.table"] = "\n".join(rows)

    # --- R3 fixed-step drift per integrator (time histories) ---
    for name in INTEGRATORS:
        tr = trajectory(out, f"rb_r3_{name}")
        m[f"r3.{name}.dE"] = max_rel_drift(tr["E_J"])
        m[f"r3.{name}.dH"] = max_vec_rel_drift(hn(tr))
        m[f"r3.{name}.dq"] = max(abs(q - 1) for q in tr["qnorm"])
    for name in ("euler", "rk4"):
        tr = trajectory(out, f"rb_r3_{name}_norenorm")
        m[f"r3.{name}_norenorm.dq"] = max(abs(q - 1) for q in tr["qnorm"])
        m[f"r3.{name}_norenorm.t_final"] = tr["t_s"][-1]
    m["r3.h"], m["r3.t_final"] = 0.02, trajectory(out, "rb_r3_rk4")["t_s"][-1]

    # --- R4 flip ---
    tr = trajectory(out, "rb_dzhanibekov")
    m["r4.sign_changes"] = sign_changes(tr["q_radps"])
    m["r4.omega2.min"], m["r4.omega2.max"] = min(tr["q_radps"]), max(tr["q_radps"])
    m["r4.dE"], m["r4.dH"] = max_rel_drift(tr["E_J"]), max_vec_rel_drift(hn(tr))
    m["r4.t_final"] = tr["t_s"][-1]
    m["r4.h"] = meta["dt_s"]

    # --- T4 phase radius on the omega = 1 rad/s axis (north), final / initial ---
    for name in ("euler", "symplectic_euler", "rk4"):
        tr = trajectory(out, f"rb_t4_{name}")
        r = [math.hypot(p, v) for p, v in zip(tr["pN_m"], tr["vN_mps"])]
        m[f"t4.{name}.radius_ratio"] = r[-1] / r[0]
    m["t4.h"], m["t4.t_final"] = 0.05, tr["t_s"][-1]

    # --- RK4 vs RKMK4 conservation at h = 0.01 (norm_conservation.csv) ---
    with open(out / "studies" / STUDY / "norm_conservation.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["case"] == "R3":
                tag = "renorm" if r["renormalize"] == "1" else "norenorm"
                for col, k in (("max_err_qnorm", "qnorm"), ("max_err_E", "E"), ("max_err_H", "H")):
                    m[f"cons.R3.{r['integrator']}.{tag}.{k}"] = float(r[col])
                m["cons.h"] = float(r["h"])

    # --- benchmarks ---
    bench = json.loads((out / "bench.json").read_text())
    ctx = bench["context"]
    m["bench.r3run.h"] = 0.01  # step used by BM_R3Run10s in bench/bench_integrators.cpp
    m["bench.cpus"], m["bench.mhz"] = ctx["num_cpus"], ctx["mhz_per_cpu"]
    for b in bench["benchmarks"]:
        kind, integ = b["name"].split("/")[0], b["label"]
        unit = {"ns": 1.0, "us": 1e3, "ms": 1e6}[b["time_unit"]]
        ns = b["cpu_time"] * unit
        if kind == "BM_Step":
            m[f"bench.step.{integ}.ns"] = ns
        else:
            m[f"bench.r3run.{integ}.us"] = ns / 1e3
    m.update(aircraft(out))
    m.update(jsbsim(out))
    m.update(control(out))
    return m


# ---------------------------------------------------------------------------------------------- control (M7)
FT = 0.3048  # reporting boundary only; the C++ runner reports SI
CONTROL_RUNS = {"c1": "f16_ap_c1_alt100", "c2": "f16_ap_c2_alt1000", "c2n": "f16_ap_c2_alt1000_nowindup",
                "c3": "f16_ap_c3_speed50", "c4": "f16_ap_c4_disturbance"}
ELEVATOR_LIMITS = {"position_deg": 25.0, "rate_dps": 60.0}  # data/f16/f16_sl.toml [actuators.elevator]


def control_meta(out: Path, tag: str) -> dict:
    return json.loads((out / "control" / CONTROL_RUNS[tag] / "metadata.json").read_text())


def _mark(ok: bool) -> str:
    return "pass" if ok else "**FAIL**"


def control(out: Path) -> dict:
    """PRD 8.7 table: targets vs C++-measured values vs the design's linear prediction (design.py output)."""
    m: dict = {}
    meta = {tag: control_meta(out, tag) for tag in CONTROL_RUNS}
    for tag, md in meta.items():
        m[f"control.{tag}.out_of_domain"] = md["out_of_domain"]
        m[f"control.{tag}.stop"] = md["stop_reason"]
        c = md["control"]
        for k, v in c["actuator_usage"].items():
            m[f"control.{tag}.use.{k}"] = v
        if c["steps"]:
            st = c["steps"][0]
            m[f"control.{tag}.overshoot"] = st["overshoot_pct"]
            m[f"control.{tag}.settling"] = st.get("settling_2pct_s", float("nan"))
            m[f"control.{tag}.settled"] = st["settling_2pct_s_reached"]
            m[f"control.{tag}.sse_ft"] = st["steady_state_error"] / FT
            m[f"control.{tag}.cross_ft"] = st["peak_cross_deviation"] / FT
            m[f"control.{tag}.t_cmd"] = st["t_cmd_s"]
            m[f"control.{tag}.t_final"] = md["t_stop_s"]
    # --- linear predictions (design.py output, step at t = 0)
    lin = json.loads((out / "control_design" / "linear_steps.json").read_text())
    for tag, key in (("c1", "C1"), ("c3", "C3")):
        lm = lin[key]["metrics"]
        m[f"control.{tag}.lin_overshoot"] = lm["overshoot_pct"]
        m[f"control.{tag}.lin_settling"] = lm["settling_2pct_s"]
        m[f"control.{tag}.lin_peak_alt_ft"] = lm["peak_altitude_ft"]
        m[f"control.{tag}.lin_peak_speed_fps"] = lm["peak_airspeed_fps"]
    # C1 steady state depends on run length (slow altitude PI tail): linear vs nonlinear |e| 55 s after the command
    t = lin["C1"]["time_s"]
    e_lin = lambda s: abs(lin["C1"]["altitude_m"][round(s / (t[1] - t[0]))] - lin["C1"]["step_si"]) / FT
    tr = trajectory_control(out, "c1")
    t0 = m["control.c1.t_cmd"]
    for horizon in (55, 300):
        m[f"control.c1.lin_e{horizon}_ft"] = e_lin(horizon)
    h0, hc = tr["h_m"][0], tr["h_cmd_m"][-1]
    idx = min(range(len(tr["t_s"])), key=lambda i: abs(tr["t_s"][i] - (t0 + 55)))
    m["control.c1.nl_e55_ft"] = abs(tr["h_m"][idx] - hc) / FT
    # --- C2 windup comparison
    for tag in ("c2", "c2n"):
        tr2 = trajectory_control(out, tag)
        m[f"control.{tag}.final_alt_error_ft"] = abs(tr2["h_m"][-1] - tr2["h_cmd_m"][-1]) / FT
        m[f"control.{tag}.theta_cmd_max_deg"] = max(tr2["theta_cmd_rad"]) / 0.017453292519943295
        m[f"control.{tag}.throttle_max"] = max(tr2["throttle"])
    # --- C4
    for d in meta["c4"]["control"]["disturbances"]:
        key = "headwind" if d["kind"] == "wind_step" else "gust"
        for loop, scale, unit in (("altitude", FT, "ft"), ("airspeed", FT, "fps")):
            x = d[loop]
            m[f"control.c4.{key}.{loop}_peak_{unit}"] = x["peak_deviation"] / scale
            m[f"control.c4.{key}.{loop}_recovery"] = x.get("recovery_2pct_s", float("nan"))
            m[f"control.c4.{key}.{loop}_recovered"] = x["recovery_2pct_s_reached"]
        m[f"control.c4.{key}.t_s"], m[f"control.c4.{key}.t_end_s"] = d["t_s"], d["t_end_s"]
    # --- C5 margins (design.py margins.json)
    margins = json.loads((out / "control_design" / "margins.json").read_text())
    rows, worst_gm, worst_pm = [], float("inf"), float("inf")
    for name, label in (("pitch_rate", "pitch rate → elevator (PI)"), ("pitch", "pitch → pitch-rate command (P)"),
                        ("altitude", "altitude → pitch command (PI)"), ("airspeed", "airspeed → throttle (PI)")):
        r = margins[name]
        lo = r["gain_reduction_margin_db"]
        worst_gm, worst_pm = min(worst_gm, r["gm_db"]), min(worst_pm, r["pm_deg"])
        rows.append(f"| {label} | {r['crossover_rad_s']:.3f} | {r['gm_db']:.2f} | {'none found' if lo is None else format(lo, '.2f')} "
                    f"| {r['pm_deg']:.2f} | {_mark(r['gm_db'] >= 6 and r['pm_deg'] >= 45 and (lo is None or lo >= 6))} |")
    m["control.c5.table"] = "\n".join(rows)
    m["control.c5.worst_gm"], m["control.c5.worst_pm"] = worst_gm, worst_pm
    # --- target table (PRD 8.7; thresholds repeated here only to label rows, asserted in test_f16_autopilot.cpp)
    g = m.__getitem__
    t = []
    t.append(f"| C1 +100 ft | overshoot | ≤ 10 % | {g('control.c1.lin_overshoot'):.2f} % | {g('control.c1.overshoot'):.2f} % | {_mark(g('control.c1.overshoot') <= 10)} |")
    t.append(f"| C1 | 2 % settling | ≤ 25 s | {g('control.c1.lin_settling'):.1f} s | {g('control.c1.settling'):.2f} s | {_mark(g('control.c1.settled') and g('control.c1.settling') <= 25)} |")
    t.append(f"| C1 | steady-state \\|e\\| (last 10 s of {g('control.c1.t_final'):.0f} s) | ≤ 1 ft | {g('control.c1.lin_e300_ft'):.2f} ft (at 300 s) | {g('control.c1.sse_ft'):.2f} ft | {_mark(g('control.c1.sse_ft') <= 1)} |")
    t.append(f"| C1 | airspeed deviation | ≤ 5 ft/s | {g('control.c1.lin_peak_speed_fps'):.2f} ft/s | {g('control.c1.cross_ft'):.2f} ft/s | {_mark(g('control.c1.cross_ft') <= 5)} |")
    t.append(f"| C2 +1,000 ft, anti-windup | overshoot | ≤ 15 % | n/a (saturates) | {g('control.c2.overshoot'):.2f} % | {_mark(g('control.c2.overshoot') <= 15)} |")
    t.append(f"| C2 | 2 % settling | ≤ 120 s | n/a | {g('control.c2.settling'):.1f} s | {_mark(g('control.c2.settled') and g('control.c2.settling') <= 120)} |")
    t.append(f"| C2 | final \\|e\\| (last 10 s of {g('control.c2.t_final'):.0f} s) | ≤ 10 ft | n/a | {g('control.c2.sse_ft'):.2f} ft | {_mark(g('control.c2.sse_ft') <= 10)} |")
    t.append(f"| C2 without anti-windup | overshoot / 2 % settling | anti-windup no worse on both | n/a | {g('control.c2n.overshoot'):.2f} % / {g('control.c2n.settling'):.1f} s | {_mark(g('control.c2.overshoot') <= g('control.c2n.overshoot') and (not g('control.c2n.settled') or g('control.c2.settling') <= g('control.c2n.settling')))} |")
    t.append(f"| C3 +50 ft/s | overshoot | ≤ 10 % | {g('control.c3.lin_overshoot'):.2f} % | {g('control.c3.overshoot'):.2f} % | {_mark(g('control.c3.overshoot') <= 10)} |")
    t.append(f"| C3 | 2 % settling | reported (should settle) | {g('control.c3.lin_settling'):.1f} s | {g('control.c3.settling'):.2f} s | {_mark(g('control.c3.settled'))} |")
    t.append(f"| C3 | altitude deviation | ≤ 20 ft | {g('control.c3.lin_peak_alt_ft'):.2f} ft | {g('control.c3.cross_ft'):.2f} ft | {_mark(g('control.c3.cross_ft') <= 20)} |")
    for key, label in (("headwind", "20 kt headwind change"), ("gust", "10 ft/s 1-cos gust")):
        rec = g(f"control.c4.{key}.altitude_recovery")
        t.append(f"| C4 {label} | peak altitude deviation / time to stay within 2 % of it | reported; returns within band | n/a | "
                 f"{g(f'control.c4.{key}.altitude_peak_ft'):.2f} ft / {rec:.1f} s | {_mark(g(f'control.c4.{key}.altitude_recovered'))} |")
    t.append(f"| C5 | worst gain margin / phase margin over four loops | ≥ 6 dB / ≥ 45° | {g('control.c5.worst_gm'):.2f} dB / {g('control.c5.worst_pm'):.2f}° | (design.py) | {_mark(g('control.c5.worst_gm') >= 6 and g('control.c5.worst_pm') >= 45)} |")
    m["control.table"] = "\n".join(t)
    # --- actuator usage (per-step peaks from the C++ runner)
    use = []
    for tag, label in (("c1", "C1"), ("c2", "C2 anti-windup"), ("c2n", "C2 no anti-windup"), ("c3", "C3"), ("c4", "C4")):
        u = {k: m[f"control.{tag}.use.{k}"] for k in ("elevator_position_peak_rad", "elevator_rate_peak_radps", "throttle_min", "throttle_max")}
        deg = 57.29577951308232
        use.append(f"| {label} | {u['elevator_position_peak_rad'] * deg:.2f} | {u['elevator_rate_peak_radps'] * deg:.2f} "
                   f"| {u['throttle_min']:.3f} … {u['throttle_max']:.3f} |")
    m["control.c1.rate_dps"] = m["control.c1.use.elevator_rate_peak_radps"] * 57.29577951308232
    m["control.usage.table"] = "\n".join(use)
    m["control.all_clean"] = not any(m[f"control.{k}.out_of_domain"] for k in CONTROL_RUNS)
    return m


def trajectory_control(out: Path, tag: str) -> dict[str, list[float]]:
    return read_csv(out / "control" / CONTROL_RUNS[tag] / "trajectory.csv")


def aircraft(out: Path) -> dict:
    """Published comparisons of C++ outputs; no dynamics, differentiation, or eigenanalysis in Python."""
    import tomllib

    refs = tomllib.loads((Path(__file__).resolve().parents[2] / "data/f16/reference_values.toml").read_text())
    m = {}
    rows = []
    for case in ("level_35", "level_30", "level_38"):
        result = json.loads((out / "aircraft" / f"trim_{case}.json").read_text())["trim"]
        for field, name, unit in (("alpha_rad", "alpha", "rad"), ("throttle", "throttle", "1"),
                                  ("elevator_deg", "elevator", "deg")):
            ours, published = result[field], refs["TR1"][case][name]["value"]
            difference = ours - published
            m[f"trim.{case}.{name}.difference"] = difference
            rows.append(f"| {case} | {name} [{unit}] | {ours:.10g} | {published:.10g} | {difference:.6g} | "
                        f"{abs(difference / published):.3%} |")
    m["trim.table"] = "\n".join(rows)
    linear = json.loads((out / "aircraft" / "linear_level_30.json").read_text())
    rows = []
    for mode in linear["modes"]:
        key = "roll_subsidence" if mode["name"] == "roll" else mode["name"]
        p = refs["LN1"]["level_30"][key]["value"]
        if not isinstance(p, list):
            p = [p, 0]
        real, imag = mode["real"], mode["imag"]
        delta = complex(real - p[0], imag - p[1])
        rel = abs(delta) / abs(complex(*p))
        rows.append(f"| {mode['name']} | {real:.9g} + {imag:.9g}i | {p[0]} + {p[1]}i | "
                    f"{delta.real:.5g} + {delta.imag:.5g}i | {rel:.3%} | {mode['wn']:.6g} | "
                    f"{mode['zeta']:.6g} | {mode['tau']:.6g} |")
    m["linear.modes.table"] = "\n".join(rows)
    rows = []
    for source in ("lateral_textbook", "longitudinal_tutorial"):
        reference = refs["LN1"]["level_30"][source]
        for i, row in enumerate(reference["A"]["value"]):
            for j, published in enumerate(row):
                matrix = "longitudinal_A" if source == "longitudinal_tutorial" else "lateral_A"
                ours = linear["sl_comparison"][matrix][i][j]
                diff = ours - published
                rel = f"{abs(diff / published):.3%}" if published else "absolute"
                rows.append(f"| {source} A[{i},{j}] | {ours:.10g} | {published:.10g} | {diff:.6g} | {rel} |")
        if source == "longitudinal_tutorial":
            for i, row in enumerate(reference["B"]["value"]):
                for j, published in enumerate(row):
                    ours = linear["sl_comparison"]["longitudinal_B"][i][j]
                    diff = ours - published
                    rel = f"{abs(diff / published):.3%}" if published else "absolute"
                    rows.append(f"| {source} B[{i},{j}] | {ours:.10g} | {published:.10g} | {diff:.6g} | {rel} |")
    m["linear.matrices.table"] = "\n".join(rows)
    m["linear.ln3.table"] = "\n".join(f"| {p['epsilon']:.1e} | {p['change']:.6g} |" for p in linear["sweep"])
    m["linear.ln3.epsilon"] = linear["chosen_epsilon"]
    m["linear.ln3.flat"] = "; ".join(f"ε={p['epsilon']:.1e}: {p['change']:.6g}" for p in linear["flat_region"])
    m["linear.ln2.altitude"] = linear["response_operating_point"]["h_m"]
    big, small = linear["responses"]
    m["linear.ln2.ratio"] = big["max_error"] / small["max_error"]
    m["linear.ln2.status"] = "pass" if all(not r["cells_changed"] and not r["out_of_domain"]
                                           and r["max_abs_aileron"] == r["max_abs_rudder"] == 0
                                           for r in linear["responses"]) and 70 <= m["linear.ln2.ratio"] <= 130 else "FAIL"
    m["linear.ln2.monitor"] = "; ".join(
        f"{r['step_deg']:.2g} deg: active_cells_changed={r['cells_changed']}, out_of_domain={r['out_of_domain']}, "
        f"max |da|={r['max_abs_aileron']:.3g} rad, max |dr|={r['max_abs_rudder']:.3g} rad"
        for r in linear["responses"])
    rows = []
    for response in linear["responses"]:
        for change in response["cell_changes"]:
            reason = change["reason"].replace("|", r"\|")
            verdict = "excluded" if change["ignored"] else "FAIL"
            rows.append(f"| {response['step_deg']:.2g} | {change['table']} | {change['from']} → {change['to']} | "
                        f"{change['t']:.3g} | {verdict}: {reason} |")
    m["linear.ln2.changes.table"] = "\n".join(rows)
    return m


SIGNAL_ORDER = ("V", "alpha", "q", "theta", "h", "beta", "p", "r", "phi")
SIGNAL_LABEL = {"V": "ΔV [ft/s]", "alpha": "Δα [deg]", "q": "Δq [deg/s]", "theta": "Δθ [deg]", "h": "Δh [ft]",
                "beta": "Δβ [deg]", "p": "Δp [deg/s]", "r": "Δr [deg/s]", "phi": "Δφ [deg]"}
CASE_LABEL = {"elev_doublet": "elevator doublet", "ail_pulse": "aileron pulse", "rud_doublet": "rudder doublet"}
BUDGET_ROWS = [  # (key in compare.json budget, label) -- order = presentation order
    ("jsbsim_dt_halving", "JSBSim step 1/2400 → 1/4800 s"),
    ("our_dt_halving", "Our RK4 step 1/2400 → 1/1200 s"),
    ("gravity_constant", "Our gravity: inverse-square → constant g0"),
    ("atmosphere", "Atmosphere ρ mismatch (sensitivity probe × measured ρ difference)"),
    ("curvature", "Earth curvature (JSBSim no-input drift removed)"),
    ("residual_trim_mismatch", "Residual trim mismatch (ours started from JSBSim's trim)"),
    ("sampling_alignment", "Sampling alignment (time mismatch × peak rate)"),
    ("table_edges", "Table edges (monitor / domain assertions)"),
]


def _pct(v):
    return "—" if v is None else f"{100 * v:.3f}"


def jsbsim(out: Path) -> dict:
    """Flatten out/jsbsim/compare.json (written by jsbsim_compare.compare) into report keys. Pure formatting."""
    c = json.loads((out / "jsbsim" / "compare.json").read_text())
    m, crit = {}, c["criteria"]
    rows = []
    for r in c["J1"]:
        j, o = r["jsbsim"], r["ours"]
        rows.append(f"| {r['h_ft']} | {j['udot_fps2']:.3e} | {j['wdot_fps2']:.3e} | {j['qdot_rps2']:.3e} | "
                    f"{r['wdot_from_rho_fps2']:.3e} | {o['wdot_fps2']:.1e} | {'pass' if r['pass'] else 'FAIL'} |")
    m["jsbsim.j1.table"] = "\n".join(rows)
    m["jsbsim.j1.status"] = "pass" if all(r["pass"] for r in c["J1"]) else "FAIL"
    rows = []
    for r in c["J2"]:
        for k, name in (("alpha_deg", "α [deg]"), ("de_deg", "δe [deg]"), ("throttle", "throttle [-]")):
            rows.append(f"| {r['h_ft']} | {name} | {r[k]['ours']:.8g} | {r[k]['jsbsim']:.8g} | {r[k]['rel']:.2e} | "
                        f"{'pass' if r['pass'] else 'FAIL'} |")
    m["jsbsim.j2.table"] = "\n".join(rows)
    m["jsbsim.j2.status"] = "pass" if all(r["pass"] for r in c["J2"]) else "FAIL"
    rows, worst, where = [], 0.0, ""
    for case, r in c["J3"].items():
        for s in SIGNAL_ORDER:
            v, peak = r["metrics"][s], r["peak"][s]
            if v is None:
                rows.append(f"| {CASE_LABEL[case]} | {SIGNAL_LABEL[s]} | {peak:.1e} | not excited | — | — |")
                continue
            self_ = c["budget"]["jsbsim_dt_halving"][case]["self"][s]
            rows.append(f"| {CASE_LABEL[case]} | {SIGNAL_LABEL[s]} | {peak:.4g} | {100 * v:.3f} | {100 * self_:.3f} | "
                        f"{'pass' if v <= crit['J3_rel'] else 'FAIL'} |")
            if v > worst:
                worst, where = v, f"{CASE_LABEL[case]}, {SIGNAL_LABEL[s]}"
    m["jsbsim.j3.table"] = "\n".join(rows)
    m["jsbsim.j3.worst"], m["jsbsim.j3.worst_where"] = 100 * worst, where
    m["jsbsim.j3.status"] = "pass" if all(r["pass"] for r in c["J3"].values()) else "FAIL"
    # budget: for each source, the largest |change of the J3 metric| over the three cases, per signal (percentage points)
    def over_cases(entries, field):
        return {s: max([abs(e[field][s]) for e in entries if e[field][s] is not None] or [None]) if any(
            e[field][s] is not None for e in entries) else None for s in SIGNAL_ORDER}

    base = {s: max([r["metrics"][s] for r in c["J3"].values() if r["metrics"][s] is not None] or [None]) if any(
        r["metrics"][s] is not None for r in c["J3"].values()) else None for s in SIGNAL_ORDER}
    rows = [f"| **Baseline J3 metric** | " + " | ".join(_pct(base[s]) for s in SIGNAL_ORDER) + " |"]
    for key, label in BUDGET_ROWS:
        eff = over_cases(list(c["budget"][key].values()), "effect")
        rows.append(f"| {label} | " + " | ".join(_pct(eff[s]) for s in SIGNAL_ORDER) + " |")
    rem = over_cases([{"e": v["metrics"]} for v in c["budget"]["all_explained"].values()], "e")
    rows.append("| **Unexplained remainder** (all of the above applied together) | " + " | ".join(_pct(rem[s]) for s in SIGNAL_ORDER) + " |")
    m["jsbsim.budget.table"] = "\n".join(rows)
    left = [(v, case, sig) for case, e in c["budget"]["all_explained"].items() for sig, v in e["metrics"].items() if v is not None]
    worst_left, worst_case, worst_sig = max(left)
    m["jsbsim.remainder.worst"] = 100 * worst_left
    m["jsbsim.remainder.where"] = f"{CASE_LABEL[worst_case]}, {SIGNAL_LABEL[worst_sig]}"
    m["jsbsim.remainder.second"] = 100 * max(v for v, case, sig in left if (case, sig) != (worst_case, worst_sig))
    m["jsbsim.remainder.steps"] = 100 * max(
        abs(e["effect"][s]) for key in ("jsbsim_dt_halving", "our_dt_halving") for e in [c["budget"][key][worst_case]]
        for s in SIGNAL_ORDER if e["effect"][s] is not None)
    a = c["atmosphere"]
    m["jsbsim.atm.table"] = "\n".join(
        f"| {name} | {a['ours'][k]:.9g} | {a['jsbsim'][k]:.9g} | {a['rel'][k]:.2e} |"
        for k, name in (("rho", "ρ [kg/m³]"), ("T", "T [K]"), ("a", "a [m/s]")))
    m["jsbsim.atm.h_ft"], m["jsbsim.atm.probe_rel_rho"] = a["h_ft"], a["probe_rel_rho"]
    cv = c["curvature"]
    rows = []
    for t, an, hold in zip(cv["t_s"], cv["analytic"], cv["jsbsim_hold"]):
        for k, name in (("h_ft", "Δh [ft]"), ("theta_deg", "Δθ [deg]"), ("V_fps", "ΔV [ft/s]")):
            rows.append(f"| {t:g} | {name} | {an[k]:.4g} | {hold[k]:.4g} |")
    m["jsbsim.curv.table"] = "\n".join(rows)
    m["jsbsim.sampling.mismatch"] = max(v["time_mismatch_s"] for v in c["budget"]["sampling_alignment"].values())
    return m
