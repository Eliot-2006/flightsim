"""JSBSim matched-model generator checks (PRD 8.6, plan Task 6.2).

The S&L evaluation in `sl_reference` is written directly from data/f16/SOURCES.md and
f16_sl.toml and deliberately shares NO code with jsbsim_compare.gen_xml: the point is two
independent encodings of the same equations, compared through JSBSim itself.
"""
import math
import tomllib
import xml.etree.ElementTree as ET
from pathlib import Path

import jsbsim
import numpy as np
import pytest

from jsbsim_compare import gen_xml

REPO = Path(__file__).resolve().parents[2]
TOML = tomllib.loads((REPO / "data/f16/f16_sl.toml").read_text())
REL = 1e-9  # task acceptance: aero point check and mass/inertia at <= 1e-9 relative


@pytest.fixture(scope="module")
def root(tmp_path_factory):
    r = tmp_path_factory.mktemp("jsb")
    gen_xml.generate(r)
    return r


def make(root, xcg=0.35):
    gen_xml.generate(root, xcg=xcg)  # overwrite aircraft file with requested CG
    fdm = jsbsim.FGFDMExec(str(root))
    fdm.set_debug_level(0)
    assert fdm.load_planet(str(root / "planet/flat_nonrot.xml"), False)
    assert fdm.load_model("f16_sl")
    fdm["simulation/gravity-model"] = 0
    return fdm


def close(a, b, rel=REL):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)


# ---- independent S&L reference (from TOML + SOURCES.md only) -------------------------------
def bilin(ax0, ax1, vals, x, y):
    ax0, ax1, vals = np.asarray(ax0, float), np.asarray(ax1, float), np.asarray(vals, float)
    i = int(np.clip(np.searchsorted(ax0, x) - 1, 0, len(ax0) - 2))
    j = int(np.clip(np.searchsorted(ax1, y) - 1, 0, len(ax1) - 2))
    tx = (x - ax0[i]) / (ax0[i + 1] - ax0[i])
    ty = (y - ax1[j]) / (ax1[j + 1] - ax1[j])
    return ((1 - tx) * (1 - ty) * vals[i, j] + tx * (1 - ty) * vals[i + 1, j]
            + (1 - tx) * ty * vals[i, j + 1] + tx * ty * vals[i + 1, j + 1])


def lin(ax, vals, x):
    return float(np.interp(x, ax, vals))  # in-domain only


def sl_reference(a_deg, b_deg, p, q, r, de, da, dr, V, qbar, xcg):
    A = TOML["aero"]
    al = A["cx"]["alpha"]
    cq, pb, rb = q * TOML["mass_geometry"]["cbar"]["value"] / (2 * V), None, None
    S, b, cbar = (TOML["mass_geometry"][k]["value"] for k in ("S", "b", "cbar"))
    xref = TOML["mass_geometry"]["xcg_ref"]["value"]
    pb, rb = p * b / (2 * V), r * b / (2 * V)
    dmp = dict(zip(A["damping"]["derivative"], [lin(al, [row[k] for row in A["damping"]["values"]], a_deg) for k in range(9)]))
    sgn = math.copysign(1.0, b_deg) if b_deg != 0 else 0.0
    cx = bilin(al, A["cx"]["elevator"], A["cx"]["values"], a_deg, de) + dmp["CXq"] * cq
    cy = (A["cy"]["beta_gain"] * b_deg + A["cy"]["aileron_gain"] * da / A["control_normalization"]["aileron"]
          + A["cy"]["rudder_gain"] * dr / A["control_normalization"]["rudder"] + dmp["CYr"] * rb + dmp["CYp"] * pb)
    cz = (lin(al, A["cz"]["values"], a_deg) * (1 - (b_deg / A["cz"]["beta_scale"]) ** 2)
          + A["cz"]["elevator_gain"] * de / A["cz"]["elevator_scale"] + dmp["CZq"] * cq)
    def eff(name, ctrl, norm):
        return bilin(al, A[name]["beta"], A[name]["values"], a_deg, b_deg) * ctrl / norm
    na, nr = A["control_normalization"]["aileron"], A["control_normalization"]["rudder"]
    cl = (sgn * bilin(al, A["cl"]["abs_beta"], A["cl"]["values"], a_deg, abs(b_deg))
          + eff("dlda", da, na) + eff("dldr", dr, nr) + dmp["Clr"] * rb + dmp["Clp"] * pb)
    cm = bilin(al, A["cm"]["elevator"], A["cm"]["values"], a_deg, de) + dmp["Cmq"] * cq + cz * (xref - xcg)
    cn = (sgn * bilin(al, A["cn"]["abs_beta"], A["cn"]["values"], a_deg, abs(b_deg))
          + eff("dnda", da, na) + eff("dndr", dr, nr) + dmp["Cnr"] * rb + dmp["Cnp"] * pb
          - cy * (xref - xcg) * cbar / b)
    qs = qbar * S
    return np.array([qs * cx, qs * cy, qs * cz, qs * b * cl, qs * cbar * cm, qs * b * cn])


def sl_thrust(throttle, h, mach):
    E = TOML["engine"]
    g = E["throttle_gearing"]
    P = g["low_slope"] * throttle if throttle <= g["breakpoint"] else g["high_slope"] * throttle + g["high_intercept"]
    t = {k: bilin(E[k]["altitude"], E[k]["mach"], E[k]["values"], h, mach) for k in ("idle", "military", "maximum")}
    if P < 50:
        return t["idle"] + (t["military"] - t["idle"]) * P / 50
    return t["military"] + (t["maximum"] - t["military"]) * (P - 50) / 50


# ---- tests ---------------------------------------------------------------------------------
STATES = [  # (V fps, h ft, alpha, beta, p, q, r deg & deg/s, de, da, dr deg, xcg)
    (502, 0, 3.0, 0.0, 0, 0, 0, 0, 0, 0, 0.35),
    (502, 0, 2.1, 0.0, 0.0, 0.0, 0.0, -0.76, 0, 0, 0.35),
    (300, 10000, 12.5, 4.0, 5.0, -3.0, 2.0, -6.0, 3.0, -4.0, 0.35),
    (450, 5000, 22.0, -7.5, -8.0, 4.0, -6.0, 9.0, -5.0, 10.0, 0.35),
    (350, 2000, 8.0, -2.0, 10.0, 6.0, 3.0, 5.0, 8.0, -8.0, 0.30),  # xcg != xcg_ref: geometric CG transfer
    (350, 2000, 18.0, 5.5, -4.0, -2.0, 4.0, -10.0, -2.0, 6.0, 0.38),
]


@pytest.mark.parametrize("st", STATES)
def test_aero_point_check(root, st):
    V, h, a, b, p, q, r, de, da, dr, xcg = st
    fdm = make(root, xcg)
    fdm["ic/vt-fps"], fdm["ic/h-sl-ft"] = V, h
    fdm["ic/alpha-deg"], fdm["ic/beta-deg"] = a, b
    fdm["ic/p-rad_sec"], fdm["ic/q-rad_sec"], fdm["ic/r-rad_sec"] = (math.radians(x) for x in (p, q, r))
    fdm["aero/de-direct-deg"], fdm["aero/da-direct-deg"], fdm["aero/dr-direct-deg"] = de, da, dr
    assert fdm.run_ic()
    got = np.array([fdm[k] for k in ("forces/fbx-aero-lbs", "forces/fby-aero-lbs", "forces/fbz-aero-lbs",
                                      "moments/l-aero-lbsft", "moments/m-aero-lbsft", "moments/n-aero-lbsft")])
    # inputs read back from JSBSim so IC round-trip rounding is not part of the check
    exp = sl_reference(fdm["aero/alpha-deg"], fdm["aero/beta-deg"], fdm["velocities/p-aero-rad_sec"],
                       fdm["velocities/q-aero-rad_sec"], fdm["velocities/r-aero-rad_sec"],
                       fdm["aero/de-total-deg"], fdm["aero/da-total-deg"], fdm["aero/dr-total-deg"],
                       fdm["velocities/vt-fps"], fdm["aero/qbar-psf"], xcg)
    assert fdm["aero/de-total-deg"] == pytest.approx(de, abs=1e-12)
    scale = np.max(np.abs(exp))
    for g, e, n in zip(got, exp, "X Y Z L M N".split()):
        # per-component rel 1e-9 with a floor of 1e-9 of the largest component (components can cancel to ~0)
        assert abs(g - e) <= REL * max(abs(e), scale * 1e-3), (n, g, e)


@pytest.mark.parametrize("throttle,h,V", [(0.2, 0, 502), (0.5, 10000, 502), (0.77, 5000, 400), (0.95, 3000, 450)])
def test_thrust_point_check(root, throttle, h, V):
    fdm = make(root)
    fdm["ic/vt-fps"], fdm["ic/h-sl-ft"] = V, h
    fdm["fcs/throttle-cmd-norm[0]"] = throttle
    assert fdm.run_ic()
    fdm["fcs/throttle-cmd-norm[0]"] = throttle
    fdm.run()
    exp = sl_thrust(throttle, fdm["position/h-sl-ft"], fdm["velocities/mach"])
    assert close(fdm["forces/fbx-external-lbs"], exp), (fdm["forces/fbx-external-lbs"], exp)


def test_mass_and_inertia(root):
    fdm = make(root)
    assert fdm.run_ic() or True
    assert close(fdm["inertia/mass-slugs"], 1 / TOML["mass_geometry"]["reciprocal_mass"]["value"])
    m = TOML["mass_geometry"]
    assert close(fdm["inertia/ixx-slugs_ft2"], m["Ixx"]["value"])
    assert close(fdm["inertia/iyy-slugs_ft2"], m["Iyy"]["value"])
    assert close(fdm["inertia/izz-slugs_ft2"], m["Izz"]["value"])
    # spike FINDINGS Q8: with negated_crossproduct_inertia="false" the property is the tensor entry -Ixz
    assert close(fdm["inertia/ixz-slugs_ft2"], -m["Ixz"]["value"])


def test_xml_tables_equal_toml(root):
    tree = ET.parse(root / "aircraft/f16_sl/f16_sl.xml")
    funcs = {f.get("name"): f for f in tree.iter("function") if f.get("name")}
    def table(name):
        t = funcs[name].find(".//table")
        rows = [list(map(float, ln.split())) for ln in t.find("tableData").text.strip().splitlines()]
        return rows
    A = TOML["aero"]
    for key, name in [("cx", "aero/table/cx"), ("cm", "aero/table/cm"), ("cl", "aero/table/cl"), ("cn", "aero/table/cn"),
                      ("dlda", "aero/table/dlda"), ("dldr", "aero/table/dldr"), ("dnda", "aero/table/dnda"), ("dndr", "aero/table/dndr")]:
        rows = table(name)
        axes = [k for k in A[key] if k in ("elevator", "abs_beta", "beta")][0]
        assert rows[0] == [float(x) for x in A[key][axes]], name
        assert [r[0] for r in rows[1:]] == [float(x) for x in A[key]["alpha"]], name
        assert [r[1:] for r in rows[1:]] == [list(map(float, v)) for v in A[key]["values"]], name
    rows = table("aero/table/cz")
    assert [r[1] for r in rows] == [float(x) for x in A["cz"]["values"]]
    for k, name in enumerate(A["damping"]["derivative"]):
        rows = table("aero/table/" + name.lower())
        assert [r[1] for r in rows] == [row[k] for row in A["damping"]["values"]], name
    for k in ("idle", "military", "maximum"):
        rows = table("thrust/table/" + k)
        assert [r[1:] for r in rows[1:]] == [list(map(float, v)) for v in TOML["engine"][k]["values"]], k


def test_committed_reference_is_current(tmp_path):
    gen_xml.generate(tmp_path)
    for rel in ("aircraft/f16_sl/f16_sl.xml", "engine/dummy.xml", "engine/direct.xml", "planet/flat_nonrot.xml"):
        assert (tmp_path / rel).read_text() == (REPO / "reference/jsbsim" / rel).read_text(), rel
