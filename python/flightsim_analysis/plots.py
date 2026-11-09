"""Plotly figures. Each entry of FIGURES is one report figure; add a figure by writing a builder and appending
a Figure(...). Builders take the `out/` directory and return a plotly Figure. Traces carry `meta` = categorical
slot index so the page can re-colour them for dark mode (see report.py)."""

import math
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .metrics import CONTROL_RUNS, ELEVATOR_LIMITS, FT, INTEGRATORS, STUDY, control_meta, hn, read_csv, trajectory, trajectory_control

# Categorical slots (dataviz palette, slots 2-4 darkened to reach 3:1 on the light surface; report.py swaps in the dark steps). Colour follows
# the integrator, never its rank in a panel.
LIGHT = ["#2a78d6", "#e0561a", "#0f8a5f", "#b8860b", "#6250d6"]  # each >= 3:1 on #fcfcfb (test_report.py)
START_MARKER = "#888888"
SLOT = {name: i for i, name in enumerate(INTEGRATORS)}
DASH = {"euler": "solid", "rk4": "solid", "symplectic_euler": "dash", "lie_euler": "dot", "rkmk4": "dashdot"}
LABEL = {"euler": "Euler", "rk4": "RK4", "symplectic_euler": "Symplectic Euler", "lie_euler": "Lie–Euler",
         "rkmk4": "RKMK4"}
NOMINAL_ORDER = {"euler": 1, "symplectic_euler": 1, "lie_euler": 1, "rk4": 4, "rkmk4": 4}
FLOOR = 1e-17  # zeros cannot be drawn on a log axis; drawn at this floor and noted in captions


@dataclass(frozen=True)
class Figure:
    id: str
    title: str
    caption: str  # what the figure demonstrates
    test_ref: str  # which test asserts it
    build: Callable[[Path], go.Figure]
    headline: bool = False  # exported to docs/img for the README


def _style(fig: go.Figure, height: int) -> go.Figure:
    fig.update_layout(template="plotly_white", height=height, margin=dict(l=70, r=20, t=40, b=60),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font=dict(family="system-ui, sans-serif", size=13),
                      legend=dict(orientation="h", y=-0.18))
    fig.update_xaxes(exponentformat="power")
    fig.update_yaxes(exponentformat="power")
    return fig


def _line(name: str, **kw) -> dict:
    return dict(meta=SLOT[name], line=dict(color=LIGHT[SLOT[name]], dash=DASH[name], width=2), **kw)


# ---------------------------------------------------------------- convergence / work-precision
# Metric plotted per case: the one the PRD 10.3 table asserts for most methods.
PRIMARY = {"T2": "err_pos", "T3": "err_pos", "T4": "err_pos", "R1": "err_att", "R2": "err_att", "R3": "err_omega"}
UNITS = {"err_pos": "position error [m]", "err_att": "attitude error [rad]", "err_omega": "rate error [rad/s]"}
OVERRIDE = {("R1", "euler"): "err_qnorm"}  # Euler on R1 is asserted on the unprojected quaternion norm


def _study(out: Path):
    return read_csv(out / "studies" / STUDY / "convergence.csv")


def _series(out: Path):
    """{(case, integrator): (xs_h, xs_nfev, ys)} for the primary metric, zeros dropped."""
    with open(out / "studies" / STUDY / "convergence.csv") as f:
        import csv

        rows = list(csv.DictReader(f))
    series: dict = {}
    for r in rows:
        key = (r["case"], r["integrator"])
        if r["case"] not in PRIMARY:
            continue
        if r["case"] == "R1" and r["integrator"] == "euler" and r["renormalize"] != "0":
            continue
        y = float(r[OVERRIDE.get(key, PRIMARY[r["case"]])])
        if y > 0:
            series.setdefault(key, ([], [], []))
            series[key][0].append(float(r["h"]))
            series[key][1].append(float(r["nfev"]))
            series[key][2].append(y)
    return series


def _small_multiples(out: Path, xi: int, xtitle: str, slopes: bool) -> go.Figure:
    series = _series(out)
    cases = list(PRIMARY)
    fig = make_subplots(rows=2, cols=3, subplot_titles=cases, horizontal_spacing=0.09, vertical_spacing=0.2)
    seen = set()
    for (case, integ), data in series.items():
        r, c = divmod(cases.index(case), 3)
        x, y = data[xi], data[2]
        metric = OVERRIDE.get((case, integ), PRIMARY[case]).replace("err_", "")
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", name=LABEL[integ], showlegend=integ not in seen,
                                 legendgroup=integ, marker=dict(size=7, color=LIGHT[SLOT[integ]]),
                                 hovertemplate=f"{LABEL[integ]} ({metric})<br>%{{x:.3g}}, %{{y:.3g}}<extra></extra>",
                                 **_line(integ)), row=r + 1, col=c + 1)
        seen.add(integ)
        if slopes:  # dashed grey reference slope through the geometric middle of the data
            p, k = NOMINAL_ORDER[integ], len(x) // 2
            xs = [x[0], x[-1]]
            ys = [y[k] * (xi_ / x[k]) ** p for xi_ in xs]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", showlegend=False, hoverinfo="skip",
                                     line=dict(color="rgba(128,128,128,0.7)", width=1, dash="dash"), meta="ref"),
                          row=r + 1, col=c + 1)
    for i, case in enumerate(cases):
        r, c = divmod(i, 3)
        fig.update_xaxes(type="log", title_text=xtitle if r == 1 else None, row=r + 1, col=c + 1)
        fig.update_yaxes(type="log", title_text=UNITS[PRIMARY[case]], row=r + 1, col=c + 1)
    return _style(fig, 620)


def convergence(out: Path) -> go.Figure:
    return _small_multiples(out, 0, "step size h [s]", slopes=True)


def work_precision(out: Path) -> go.Figure:
    return _small_multiples(out, 1, "derivative evaluations [-]", slopes=False)


# ---------------------------------------------------------------- time histories
def drift(out: Path) -> go.Figure:
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Energy", "Angular momentum |H_N|"),
                        horizontal_spacing=0.12)
    for name in INTEGRATORS:
        tr = trajectory(out, f"rb_r3_{name}")
        e0, h0 = tr["E_J"][0], hn(tr)[0]
        de = [max(abs(e / e0 - 1), FLOOR) for e in tr["E_J"]]
        dh = [max(math.dist(h, h0) / math.hypot(*h0), FLOOR) for h in hn(tr)]
        for col, y in ((1, de), (2, dh)):
            fig.add_trace(go.Scatter(x=tr["t_s"], y=y, name=LABEL[name], legendgroup=name, showlegend=col == 1,
                                     mode="lines", **_line(name)), row=1, col=col)
    fig.update_xaxes(title_text="time t [s]")
    fig.update_yaxes(type="log", title_text="relative drift [-]")
    return _style(fig, 430)


def qnorm(out: Path) -> go.Figure:
    fig = go.Figure()
    runs = [("euler", "rb_r3_euler_norenorm", "Euler, no renormalization"),
            ("rk4", "rb_r3_rk4_norenorm", "RK4, no renormalization"),
            ("rkmk4", "rb_r3_rkmk4", "RKMK4 (no renormalization needed)"),
            ("euler", "rb_r3_euler", "Euler, renormalized"),
            ("rk4", "rb_r3_rk4", "RK4, renormalized")]
    for integ, run, label in runs:
        tr = trajectory(out, run)
        y = [max(abs(q - 1), FLOOR) for q in tr["qnorm"]]
        style = _line(integ)
        if "renormalized" in label:
            style["line"]["dash"] = "dot"
        fig.add_trace(go.Scatter(x=tr["t_s"], y=y, name=label, mode="lines", **style))
    fig.update_xaxes(title_text="time t [s]")
    fig.update_yaxes(type="log", title_text="|‖q‖ − 1| [-]")
    return _style(fig, 430)


def r4(out: Path) -> go.Figure:
    tr = trajectory(out, "rb_dzhanibekov")
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.1,
                        subplot_titles=("Body rates ω_B", "Inertial angular momentum H_N"))
    for i, (k, lab) in enumerate((("p_radps", "p"), ("q_radps", "q (intermediate axis)"), ("r_radps", "r"))):
        fig.add_trace(go.Scatter(x=tr["t_s"], y=tr[k], name=lab, mode="lines",
                                 line=dict(color=LIGHT[i], width=2), meta=i), row=1, col=1)
    for i, (k, lab) in enumerate((("HN_x_kgm2ps", "H_N,x"), ("HN_y_kgm2ps", "H_N,y"), ("HN_z_kgm2ps", "H_N,z"))):
        fig.add_trace(go.Scatter(x=tr["t_s"], y=tr[k], name=lab, mode="lines",
                                 line=dict(color=LIGHT[i], width=2, dash="dash"), meta=i), row=2, col=1)
    fig.update_xaxes(title_text="time t [s]", row=2, col=1)
    fig.update_yaxes(title_text="ω [rad/s]", row=1, col=1)
    fig.update_yaxes(title_text="H_N [kg·m²/s]", row=2, col=1)
    return _style(fig, 560)


def t4(out: Path) -> go.Figure:
    fig = go.Figure()
    for name in ("euler", "symplectic_euler", "rk4"):
        tr = trajectory(out, f"rb_t4_{name}")
        fig.add_trace(go.Scatter(x=tr["pN_m"], y=tr["vN_mps"], name=LABEL[name], mode="lines", **_line(name)))
    fig.add_trace(go.Scatter(x=[1], y=[0], mode="markers", name="start", marker=dict(size=9, color=START_MARKER),
                             meta="ref"))
    fig.update_xaxes(title_text="position p_N [m]")
    fig.update_yaxes(title_text="velocity v_N [m/s]", scaleanchor="x", scaleratio=1)
    return _style(fig, 520)


def ln3(out: Path) -> go.Figure:
    data = json.loads((out / "aircraft/linear_level_30.json").read_text())
    fig = go.Figure(go.Scatter(x=[p["epsilon"] for p in data["sweep"]],
                               y=[p["change"] for p in data["sweep"]], mode="lines+markers",
                               name="max relative entry change vs 10ε", meta=0,
                               line=dict(color=LIGHT[0], width=2)))
    fig.add_vline(x=data["chosen_epsilon"], line_dash="dash", annotation_text="chosen ε = 6e-6")
    fig.update_xaxes(type="log", title_text="finite-difference ε [-]")
    fig.update_yaxes(type="log", title_text="max entry change / max matrix entry [-]")
    return _style(fig, 430)


def ln2(out: Path) -> go.Figure:
    data = json.loads((out / "aircraft/linear_level_30.json").read_text())
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Airspeed perturbation", "Angle-of-attack perturbation"))
    for slot, response in enumerate(data["responses"]):
        points = response["points"]
        for key, dash in (("nonlinear", "solid"), ("linear", "dash")):
            for col, state in ((1, 0), (2, 1)):
                fig.add_trace(go.Scatter(x=[p["t"] for p in points], y=[p[key][state] for p in points],
                                         name=f"{response['step_deg']:.2g}° {key}", showlegend=col == 1,
                                         mode="lines", meta=slot, line=dict(color=LIGHT[slot], dash=dash)),
                              row=1, col=col)
    fig.update_xaxes(title_text="time t [s]")
    fig.update_yaxes(title_text="ΔV [m/s]", row=1, col=1)
    fig.update_yaxes(title_text="Δα [rad]", row=1, col=2)
    return _style(fig, 430)


FIGURES = [
    Figure("convergence", "Convergence: error against step size",
           "Final-time error against step size h on log–log axes, one panel per validation case; dashed grey lines "
           "are reference slopes of the method's nominal order (1 or 4). Parallel lines confirm the order. Methods with the same error on a case overlap (R3: Euler with Lie–Euler, RK4 with RKMK4).",
           "validation_tests: Translational.T2Rk4Order4, T2EulerOrder1, T3Rk4Order4, T3EulerOrder1; Rotational.R1Rk4Order4, R2Rk4Order4, R2Rkmk4Order4, R3Rk4Order4, R3Rkmk4Order4 and the first-order R1–R3 tests (observed order in [0.95, 1.05] / [3.9, 4.1]); Convergence.T4OrdersForApplicableIntegrators",
           convergence, headline=True),
    Figure("work-precision", "Work–precision: error against derivative evaluations",
           "The same errors against the number of derivative evaluations. RK4 and RKMK4 cost four evaluations per "
           "step yet reach small errors with far less work than the first-order methods.",
           "No separate assertion: the same studies as the convergence figure, with the nfev column (the order tests above gate the errors; PRD §6.3)", work_precision),
    Figure("drift", "Energy and angular-momentum drift, R3 at fixed h",
           "Relative drift of rotational energy and |H_N| for the torque-free asymmetric body R3 at fixed h = 0.02 s. "
           "Zeros are drawn at 1e-17. Euler and Lie–Euler drift secularly; RK4 and RKMK4 drift by many orders of magnitude less.",
           "validation_tests: Rotational.HNConservedRk4, Rotational.R3Rkmk4Order4, Rotational.R3LieEulerOrder1 (H-error orders)",
           drift, headline=True),
    Figure("qnorm", "Quaternion-norm drift with and without renormalization",
           "|‖q‖ − 1| over time on R3. Without projection Euler and RK4 drift away from unit norm; renormalization "
           "pins it at roundoff; RKMK4 preserves the norm by construction.",
           "validation_tests: Rotational.NormDriftWithoutRenormalization; unit_tests: LieIntegrators.LieEulerPreservesNorm", qnorm),
    Figure("r4", "R4: intermediate-axis (Dzhanibekov) flip",
           "Spin about the intermediate axis with a 1e-3 rad/s perturbation: ω₂ changes sign while H_N stays "
           "constant in the inertial frame.",
           "validation_tests: Rotational.R4FlipOccursAndConserves, SpinStability.L1SpinStability", r4, headline=True),
    Figure("t4", "T4: phase portrait of the harmonic oscillator (north axis, ω = 1 rad/s)",
           "Position–velocity orbit at a deliberately coarse h = 0.05 s. Explicit Euler spirals outward, symplectic "
           "Euler stays on a bounded closed curve, RK4 stays on a closed circle.",
           "validation_tests: Translational.T4EulerEnergyGrows, Translational.T4Rk4EnergyDriftSmall; unit_tests: LieIntegrators.SymplecticEulerOscillatorEnergyBounded",
           t4),
]

FIGURES.extend([
    Figure("ln3", "LN3: finite-difference step sweep",
           "Central differences from ε=1e-2 to 1e-10: truncation, a flat region, then roundoff. Changes are maximum absolute entry differences divided by the maximum absolute entry of each matrix, taking the larger of A and B. C++ computes every Jacobian.",
           "validation_tests: F16Linear.LN3StepSweep (chosen ε vs 10ε and ε/10 ≤1e-5)", ln3),
    Figure("ln2", "LN2: linear and nonlinear elevator responses",
           "Trailing-edge-up steps of 0.5° and 0.05° over 5 s from a fresh 100 ft trim, computed in C++. The ratio and lookup-monitor verdict are reported in docs/results.md; the guard excludes only zero-multiplier lateral controls and the asserted smooth beta=0 fold.",
           "validation_tests: F16Linear.LN2QuadraticScaling (ratio [70,130] AND no active table kink crossed)", ln2),
])


# ---------------------------------------------------------------- JSBSim cross-validation (J3, discrepancy budget)
def _compare(out: Path) -> dict:
    return json.loads((out / "jsbsim" / "compare.json").read_text())


def _j3_overlay(case: str, signals: list[str]):
    def build(out: Path) -> go.Figure:
        from .metrics import SIGNAL_LABEL

        c = _compare(out)
        s = c["series"][case]
        fig = make_subplots(rows=2, cols=2, shared_xaxes=True, vertical_spacing=0.14, horizontal_spacing=0.12,
                            subplot_titles=[SIGNAL_LABEL[k] for k in signals])
        for i, k in enumerate(signals):
            r, col = divmod(i, 2)
            fig.add_trace(go.Scatter(x=s["t_s"], y=s[k]["ours"], mode="lines", name="flightsim (RK4, 1/2400 s)", legendgroup="o",
                                     showlegend=i == 0, meta=0, line=dict(color=LIGHT[0], width=3)), row=r + 1, col=col + 1)
            fig.add_trace(go.Scatter(x=s["t_s"], y=s[k]["jsbsim"], mode="lines", name="JSBSim (AB3, 1/2400 s)", legendgroup="j",
                                     showlegend=i == 0, meta=1, line=dict(color=LIGHT[1], width=2, dash="dash")), row=r + 1, col=col + 1)
            fig.update_yaxes(title_text=SIGNAL_LABEL[k].split(" ", 1)[1].strip("[]"), row=r + 1, col=col + 1)
            if r == 1:
                fig.update_xaxes(title_text="time t [s]", row=r + 1, col=col + 1)
        return _style(fig, 560)

    return build


def j3_budget(out: Path) -> go.Figure:
    """Worst |change of any J3 metric| per budget source, against the 2 % criterion (log axis)."""
    from .metrics import BUDGET_ROWS

    c = _compare(out)
    names, vals = [], []
    for key, label in BUDGET_ROWS:
        eff = [abs(v) for e in c["budget"][key].values() for v in e["effect"].values() if v is not None]
        names.append(label)
        vals.append(max(eff) * 100)
    rem = [v for e in c["budget"]["all_explained"].values() for v in e["metrics"].values() if v is not None]
    names.append("Unexplained remainder")
    vals.append(max(rem) * 100)
    base = [v for e in c["J3"].values() for v in e["metrics"].values() if v is not None]
    names.append("Baseline J3 difference")
    vals.append(max(base) * 100)
    vals = [max(v, 1e-9) for v in vals]  # zero effects are drawn at the floor
    fig = go.Figure(go.Bar(y=names[::-1], x=vals[::-1], orientation="h", meta=0, marker=dict(color=LIGHT[0]),
                           hovertemplate="%{y}<br>%{x:.3g} % of peak<extra></extra>"))
    fig.add_vline(x=c["criteria"]["J3_rel"] * 100, line_dash="dash", line_color=START_MARKER,
                  annotation_text="J3 criterion 2 %")
    fig.update_xaxes(type="log", range=[-9.3, 1], title_text="largest change of any J3 metric [% of the signal's peak deviation]")
    fig.update_yaxes(automargin=True)
    return _style(fig, 520).update_layout(margin=dict(l=20, r=20, t=40, b=60))


FIGURES.extend([
    Figure("j3-elevator", "J3: elevator doublet, flightsim vs JSBSim",
           "Deviation from each simulator's own trim (502 ft/s, 10,000 ft) for ±1° elevator pulses at 1–2 s and 2–3 s. "
           "The curves overlap to within the line width; the numerical differences are in docs/results.md.",
           "python/tests/test_jsbsim_compare.py and `python -m jsbsim_compare.compare --check` (max difference ≤ 2 % of peak)",
           _j3_overlay("elev_doublet", ["V", "alpha", "q", "theta"])),
    Figure("j3-aileron", "J3: aileron pulse, flightsim vs JSBSim",
           "Deviation from own trim for a 2° aileron pulse at 1–2 s: sideslip, roll rate, yaw rate and bank angle.",
           "python/tests/test_jsbsim_compare.py and `python -m jsbsim_compare.compare --check`",
           _j3_overlay("ail_pulse", ["beta", "p", "r", "phi"])),
    Figure("j3-rudder", "J3: rudder doublet, flightsim vs JSBSim",
           "Deviation from own trim for ±2° rudder pulses at 1–2 s and 2–3 s.",
           "python/tests/test_jsbsim_compare.py and `python -m jsbsim_compare.compare --check`",
           _j3_overlay("rud_doublet", ["beta", "p", "r", "phi"])),
    Figure("j3-budget", "J3 discrepancy budget",
           "Each candidate source of disagreement is toggled or bounded on its own and the change in the largest J3 "
           "difference is plotted (log axis; exact zeros are drawn at the floor). The unexplained remainder applies all "
           "explained corrections at once.",
           "python/tests/test_jsbsim_compare.py (budget structure) and `python -m jsbsim_compare.compare --check`", j3_budget),
])


# ---------------------------------------------------------------- control (M7): C1-C4, windup, actuator usage
# Slots: 0 nonlinear, 1 linear prediction / second run, 2 third run, 3 fourth run; command reference is grey.
DEG = 0.017453292519943295


def _ln(name, x, y, slot, dash="solid", **kw):
    return go.Scatter(x=x, y=y, name=name, mode="lines", meta=slot,
                      line=dict(color=LIGHT[slot], dash=dash, width=2), **kw)


def _ref(name, x, y):
    return go.Scatter(x=x, y=y, name=name, mode="lines", meta="ref", line=dict(color=START_MARKER, dash="dot", width=1.5))


def _lin_steps(out: Path) -> dict:
    return json.loads((out / "control_design" / "linear_steps.json").read_text())


def _overlay_builder(tag: str, key: str, scale: float):
    def build(out: Path) -> go.Figure:
        tr = trajectory_control(out, tag)
        lin = _lin_steps(out)[key]
        t_cmd = control_meta(out, tag)["control"]["steps"][0]["t_cmd_s"]
        t0, h0, v0 = tr["t_s"], tr["h_m"][0], tr["V_mps"][0]
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                            subplot_titles=("Altitude deviation from trim", "Airspeed deviation from trim"))
        t_lin = [t_cmd + t for t in lin["time_s"]]
        for row, (nl, ln, cmd_key, ref0) in enumerate(((tr["h_m"], lin["altitude_m"], "h_cmd_m", h0),
                                                       (tr["V_mps"], lin["airspeed_mps"], "V_cmd_mps", v0)), start=1):
            fig.add_trace(_ref("command", t0, [(c - ref0) / FT for c in tr[cmd_key]]).update(showlegend=row == 1), row=row, col=1)
            fig.add_trace(_ln("nonlinear (C++ runner)", t0, [(y - ref0) / FT for y in nl], 0, showlegend=row == 1), row=row, col=1)
            fig.add_trace(_ln("linear prediction" + (f" × {scale:g}" if scale != 1 else ""), t_lin,
                              [scale * y / FT for y in ln], 1, "dash", showlegend=row == 1), row=row, col=1)
        fig.update_xaxes(title_text="time t [s]", range=[0, min(tr["t_s"][-1], t_cmd + 100)], row=2, col=1)
        fig.update_yaxes(title_text="Δh [ft]", row=1, col=1)
        fig.update_yaxes(title_text="ΔV [ft/s]", row=2, col=1)
        for i in fig.layout.annotations:
            i.font.size = 13
        return _style(fig, 620)

    return build


def c2_windup(out: Path) -> go.Figure:
    # The windup difference lives in the first ~100 s (overshoot and the slow return); a 300 s axis hides it, so the
    # figure shows 0-120 s and adds the altitude error, where overshoot and settling are read directly.
    fig = make_subplots(rows=4, cols=1, shared_xaxes=True, vertical_spacing=0.06,
                        subplot_titles=("Altitude deviation from trim", "Altitude error h − h_cmd",
                                        "Pitch command θ_cmd", "Throttle command"))
    runs = (("c2", "anti-windup on (back-calculation)", 0, "solid"), ("c2n", "anti-windup off", 1, "dash"))
    ref_done = False
    for tag, label, slot, dash in runs:
        tr = trajectory_control(out, tag)
        h0 = tr["h_m"][0]
        if not ref_done:
            fig.add_trace(_ref("command", tr["t_s"], [(c - h0) / FT for c in tr["h_cmd_m"]]), row=1, col=1)
            ref_done = True
        fig.add_trace(_ln(label, tr["t_s"], [(h - h0) / FT for h in tr["h_m"]], slot, dash), row=1, col=1)
        err = [(h - c) / FT for h, c in zip(tr["h_m"], tr["h_cmd_m"])]
        fig.add_trace(_ln(label, tr["t_s"], err, slot, dash, showlegend=False), row=2, col=1)
        fig.add_trace(_ln(label, tr["t_s"], [t / DEG for t in tr["theta_cmd_rad"]], slot, dash, showlegend=False), row=3, col=1)
        fig.add_trace(_ln(label, tr["t_s"], tr["throttle"], slot, dash, showlegend=False), row=4, col=1)
    fig.update_xaxes(range=[0, 120])
    fig.update_xaxes(title_text="time t [s]", row=4, col=1)
    fig.update_yaxes(title_text="Δh [ft]", row=1, col=1)
    fig.update_yaxes(title_text="error [ft]", range=[-60, 60], row=2, col=1)
    fig.update_yaxes(title_text="θ_cmd [deg]", row=3, col=1)
    fig.update_yaxes(title_text="throttle [-]", row=4, col=1)
    return _style(fig, 900)


def c4_disturbance(out: Path) -> go.Figure:
    tr = trajectory_control(out, "c4")
    dist = {d["kind"]: d for d in control_meta(out, "c4")["control"]["disturbances"]}
    fig = make_subplots(rows=2, cols=2, shared_xaxes=False, vertical_spacing=0.12, horizontal_spacing=0.08,
                        subplot_titles=("20 kt headwind change: altitude", "1-cos gust: altitude",
                                        "20 kt headwind change: airspeed", "1-cos gust: airspeed"))
    for col, kind, span in ((1, "wind_step", 450), (2, "gust", 200)):
        t0 = dist[kind]["t_s"]
        idx = [i for i, t in enumerate(tr["t_s"]) if t0 - 5 <= t <= t0 + span]
        rel = [tr["t_s"][i] - t0 for i in idx]
        for row, (sig, cmd, label) in enumerate((("h_m", "h_cmd_m", "altitude"), ("V_mps", "V_cmd_mps", "airspeed")), start=1):
            dev = [(tr[sig][i] - tr[cmd][i]) / FT for i in idx]
            peak = max(abs(d) for d in dev)
            fig.add_trace(_ln(label, rel, dev, 0, showlegend=False), row=row, col=col)
            for sgn in (1, -1):  # the 2 % band the recovery time is measured against
                fig.add_trace(_ref("±2 % of peak", rel, [sgn * 0.02 * peak] * len(rel)).update(showlegend=(sgn == 1 and row == 1 and col == 1)), row=row, col=col)
            fig.update_xaxes(title_text="time since disturbance [s]" if row == 2 else None, row=row, col=col)
    fig.update_yaxes(title_text="Δh from command [ft]", row=1, col=1)
    fig.update_yaxes(title_text="ΔV from command [ft/s]", row=2, col=1)
    return _style(fig, 640)


def actuator_usage(out: Path) -> go.Figure:
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07,
                        subplot_titles=("Elevator position", "Elevator rate", "Throttle command"))
    cases = (("c1", "C1 +100 ft", 0), ("c2", "C2 +1,000 ft", 1), ("c3", "C3 +50 ft/s", 2))
    for tag, label, slot in cases:
        tr = trajectory_control(out, tag)
        t_cmd = control_meta(out, tag)["control"]["steps"][0]["t_cmd_s"]
        idx = [i for i, t in enumerate(tr["t_s"]) if t_cmd - 2 <= t <= t_cmd + 60]
        rel = [tr["t_s"][i] - t_cmd for i in idx]
        de = [tr["de_state_rad"][i] / DEG for i in idx]
        # rate = first difference of the recorded actuator state (post-processing; the C++ per-step peak is in the table)
        rate = [0.0] + [(de[i] - de[i - 1]) / (rel[i] - rel[i - 1]) for i in range(1, len(de))]
        fig.add_trace(_ln(label, rel, de, slot), row=1, col=1)
        fig.add_trace(_ln(label, rel, rate, slot, showlegend=False), row=2, col=1)
        fig.add_trace(_ln(label, rel, [tr["throttle"][i] for i in idx], slot, showlegend=False), row=3, col=1)
    span = [-2, 60]
    for row, level in ((1, ELEVATOR_LIMITS["position_deg"]), (2, ELEVATOR_LIMITS["rate_dps"])):
        for sgn in (1, -1):
            fig.add_trace(_ref("actuator limit", span, [sgn * level] * 2).update(showlegend=(row == 1 and sgn == 1)), row=row, col=1)
    for level in (0, 1):
        fig.add_trace(_ref("throttle limit", span, [level] * 2).update(showlegend=False), row=3, col=1)
    fig.update_xaxes(title_text="time since command [s]", row=3, col=1)
    fig.update_yaxes(title_text="δe [deg]", row=1, col=1)
    fig.update_yaxes(title_text="dδe/dt [deg/s]", row=2, col=1)
    fig.update_yaxes(title_text="throttle [-]", row=3, col=1)
    return _style(fig, 760)


FIGURES.extend([
    Figure("c1-overlay", "C1: +100 ft altitude step, nonlinear vs linear prediction",
           "Closed-loop response of the C++ runner (50 Hz autopilot, actuator and engine lag) against the design.py linear "
           "prediction (linear_steps.json perturbations from trim, shifted to the command time). The two agree to a fraction of a percent of the step.",
           "validation_tests: F16Autopilot.C1AltitudeStep100ft (overshoot ≤ 10 %, 2 % settling ≤ 25 s, steady-state ≤ 1 ft, airspeed deviation ≤ 5 ft/s)",
           _overlay_builder("c1", "C1", 1.0), headline=True),
    Figure("c2-overlay", "C2: +1,000 ft altitude step, nonlinear vs scaled linear prediction",
           "The saturating +1,000 ft step against the linear C1 prediction multiplied by 10 (the linear model is linear but has no saturation). "
           "The gap is the effect of the ±15° pitch-command limit, throttle limit and back-calculation.",
           "validation_tests: F16Autopilot.C2AltitudeStep1000ftAntiwindupBeatsNoAntiwindup",
           _overlay_builder("c2", "C1", 10.0)),
    Figure("c3-overlay", "C3: +50 ft/s airspeed step, nonlinear vs linear prediction",
           "The +50 ft/s step needs more than the whole throttle range at the first sample (kp·ΔV exceeds 1), so the linear prediction, which has "
           "no throttle limit, does not apply; the nonlinear run saturates the throttle and recovers slowly.",
           "validation_tests: F16Autopilot.C3AirspeedStep50fps (overshoot ≤ 10 %, altitude deviation ≤ 20 ft)",
           _overlay_builder("c3", "C3", 1.0)),
    Figure("c2-windup", "C2: anti-windup against no anti-windup",
           "The +1,000 ft step with back-calculation anti-windup (solid) and without (dashed): altitude, pitch command and throttle.",
           "validation_tests: F16Autopilot.C2AltitudeStep1000ftAntiwindupBeatsNoAntiwindup (overshoot with ≤ 15 % and strictly less than without)",
           c2_windup, headline=True),
    Figure("c4-disturbance", "C4: headwind change and vertical gust rejection",
           "Deviation from the altitude and airspeed commands after a 20 kt headwind step (left) and a 10 ft/s 1-cos updraft (right). "
           "Dotted lines are the ±2 % of peak band the recovery time is measured against.",
           "validation_tests: F16Autopilot.C4GustAndHeadwindChangeRecoverWithinTwoPercentBand",
           c4_disturbance),
    Figure("control-actuators", "Control: actuator usage",
           "Elevator position and rate and the throttle command for C1–C3, against the elevator ±25° / 60°/s and throttle [0, 1] limits.",
           "unit_tests: AircraftRunnerAutopilot (rates/limits are applied by the model); the per-step peaks are in docs/results.md",
           actuator_usage),
])
