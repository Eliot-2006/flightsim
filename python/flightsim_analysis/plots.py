"""Plotly figures. Each entry of FIGURES is one report figure; add a figure by writing a builder and appending
a Figure(...). Builders take the `out/` directory and return a plotly Figure. Traces carry `meta` = categorical
slot index so the page can re-colour them for dark mode (see report.py)."""

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .metrics import INTEGRATORS, STUDY, hn, read_csv, trajectory

# Categorical slots (dataviz reference palette, light mode; report.py swaps in the dark steps). Colour follows
# the integrator, never its rank in a panel.
LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6"]
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
    fig.add_trace(go.Scatter(x=[1], y=[0], mode="markers", name="start", marker=dict(size=9, color="#888"),
                             meta="ref"))
    fig.update_xaxes(title_text="position p_N [m]")
    fig.update_yaxes(title_text="velocity v_N [m/s]", scaleanchor="x", scaleratio=1)
    return _style(fig, 520)


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
