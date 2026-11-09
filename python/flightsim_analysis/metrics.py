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
    return m


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
