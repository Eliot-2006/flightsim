# flightsim

flightsim is a C++20 six-degree-of-freedom rigid-body simulator and an F-16 flight-dynamics model, built so that every physical claim can be traced to a test or a reproducible measurement in this repository. The generic engine (quaternion attitude, five integrators) is checked against closed-form solutions, among them torque-free asymmetric rotation through Jacobi elliptic functions, with measured convergence orders and conservation laws asserted in CI. The aircraft model, NASA TP-1538 data as tabulated by Stevens, Lewis & Johnson, is checked against published check cases, work–energy and angular-impulse balances, and an independent simulator (JSBSim) running a generated copy of the same tables. Trim, linearization and an altitude and airspeed autopilot sit on top of it. Where a reference is weaker than it looks (a single secondary source, a different model), this README says so.

> Generated from `docs/README.template.md` by `python -m flightsim_analysis.report`; every number comes from `out/metrics.json`, written by the same demo run (`scripts/demo.sh`).

## Headline numbers

The last column names the test or scenario that produces each number. Bounds are the ones a test enforces; "none" means the value is reported, not asserted. [`docs/results.md`](docs/results.md) has every table and a [traceability table](docs/results.md#traceability-readme-claims) from each claim to its test, scenario and figure.

| Result | Measured | Bound | Produced by |
|---|---|---|---|
| Convergence orders inside their bounds | {{orders.in_bounds}} of {{orders.asserted}}: Euler, Lie–Euler, symplectic Euler {{orders.first.min:.3f}}–{{orders.first.max:.3f}}; RK4, RKMK4 {{orders.fourth.min:.3f}}–{{orders.fourth.max:.3f}} | [{{orders.first.lo}}, {{orders.first.hi}}] and [{{orders.fourth.lo}}, {{orders.fourth.hi}}] | [`core_convergence.toml`](scenarios/studies/core_convergence.toml); [`test_translational.cpp`](tests/validation/test_translational.cpp), [`test_rotational.cpp`](tests/validation/test_rotational.cpp) |
| Energy drift on the torque-free asymmetric body R3 over {{r3.t_final:.0f}} s | RK4 {{r3.rk4.dE:.1e}}, RKMK4 {{r3.rkmk4.dE:.1e}}, Euler {{r3.euler.dE:.1e}} | none | [`rb_r3_*.toml`](scenarios/rb_r3_rk4.toml) |
| Quaternion-norm drift on R3 without renormalization | RK4 {{r3.rk4_norenorm.dq:.1e}}; RKMK4 {{r3.rkmk4.dq:.1e}} with no projection | none | [`rb_r3_rk4_norenorm.toml`](scenarios/rb_r3_rk4_norenorm.toml), `Rotational.NormDriftWithoutRenormalization` |
| Intermediate-axis flip R4: sign changes of the body rate q in {{r4.t_final:.0f}} s | {{r4.sign_changes}}, with H_N drift {{r4.dH:.1e}} and energy drift {{r4.dE:.1e}} | flip must occur; conservation asserted | [`rb_dzhanibekov.toml`](scenarios/rb_dzhanibekov.toml), `Rotational.R4FlipOccursAndConserves` |
| F-16 state-derivative check case A1 (published S&L values) | largest relative error {{checks.a1.max_rel:.1e}} | 1e-3 | `F16Model.A1CheckCase` in [`test_f16_model.cpp`](tests/unit/test_f16_model.cpp) |
| Work–energy balance B1, angular-impulse balance B2 (20 s, dt = 0.01 s) | {{checks.b1.residual:.1e}}, {{checks.b2.residual:.1e}} | 1e-6 | [`test_f16_consistency.cpp`](tests/validation/test_f16_consistency.cpp) |
| Level-flight trim TR1 and linear-model eigenvalues LN1 against published values | largest relative difference {{trim.max_rel_pct:.3f}} % and {{linear.modes.max_rel_pct:.2f}} % | 0.5 % (or the last printed digit) and 1 % | [`f16_trim_502_*.toml`](scenarios/f16_trim_502_xcg035.toml), [`test_f16_trim.cpp`](tests/validation/test_f16_trim.cpp), [`test_f16_linear.cpp`](tests/validation/test_f16_linear.cpp) |
| Linear against nonlinear elevator step LN2: error ratio for a step 10× smaller | {{linear.ln2.ratio:.1f}} | 70 to 130 (quadratic scaling gives 100) | `F16Linear.LN2QuadraticScaling` |
| JSBSim J3: open-loop responses, worst difference over every signal and case | {{jsbsim.j3.worst:.2f}} % of the signal's peak deviation | 2 % | [`compare.py --check`](python/jsbsim_compare/compare.py), [`scenarios/jsbsim/`](scenarios/jsbsim) |
| Autopilot C1, +100 ft step: overshoot and 2 % settling time | {{control.c1.overshoot:.2f}} %, {{control.c1.settling:.1f}} s | 10 %, 25 s | [`f16_ap_c1_alt100.toml`](scenarios/f16_ap_c1_alt100.toml), `F16Autopilot.C1AltitudeStep100ft` |
| Autopilot C2, +1,000 ft step: settling time with and without anti-windup | {{control.c2.settling:.1f}} s and {{control.c2n.settling:.1f}} s | 120 s, and no worse than without | [`f16_ap_c2_alt1000.toml`](scenarios/f16_ap_c2_alt1000.toml), `F16Autopilot.C2AltitudeStep1000ftAntiwindupBeatsNoAntiwindup` |
| Stock JSBSim `f16` against ours at the same trim, thrust (sea level; 10,000 ft) | {{j4.h0.thrust_diff_pct:+.0f}} %; {{j4.h10000.thrust_diff_pct:+.0f}} % (a model-form difference) | none | [`stock_f16.py`](python/jsbsim_compare/stock_f16.py) |
| F-16 RK4 step; 60 s closed-loop run | {{bench.f16.BM_F16Rk4Step.ns:.0f}} ns; {{bench.f16.realtime_factor:.0f}}× real time | PRD §11 targets, reported not gated | [`bench_f16.cpp`](bench/bench_f16.cpp) |

## Four figures

### Five integrators converge at their expected orders

![Log–log error against step size for T2, T3, T4, R1, R2 and R3, with reference slopes of 1 and 4](docs/img/convergence.png)

Error at the final time against step size for the closed-form cases. Euler, Lie–Euler and symplectic Euler follow the slope-1 reference; RK4 and RKMK4 follow slope 4 until round-off. The fitted orders are asserted to lie in [{{orders.first.lo}}, {{orders.first.hi}}] and [{{orders.fourth.lo}}, {{orders.fourth.hi}}] (PRD §10.3). Produced by [`core_convergence.toml`](scenarios/studies/core_convergence.toml); the interpolation-order caveat for the aircraft is in [`docs/results.md`](docs/results.md).

### The intermediate-axis flip conserves angular momentum

![Body rates of a rigid body spinning about its intermediate axis, flipping sign repeatedly while the inertial angular momentum stays constant](docs/img/r4.png)

Spin about the intermediate axis of an asymmetric body is unstable, so the body rate q flips sign {{r4.sign_changes}} times in {{r4.t_final:.0f}} s while the angular momentum in the inertial frame, H_N, stays constant to {{r4.dH:.1e}} (RK4, h = {{r4.h}} s). The linearized spin eigenvalues explain the instability and match the closed form ([`SpinStability.L1SpinStability`](tests/validation/test_spin_stability.cpp)). Produced by [`rb_dzhanibekov.toml`](scenarios/rb_dzhanibekov.toml).

### Our model against JSBSim on the same data

![Elevator doublet response, flightsim and JSBSim overlaid for airspeed, angle of attack, pitch rate and pitch angle](docs/img/j3-elevator.png)

JSBSim flies an aircraft generated from the same data file as ours, on a non-rotating planet with the same gravity and atmosphere, and both start from their own trim. The overlay is the elevator doublet at 10,000 ft; the worst difference over all three maneuvers and nine signals is {{jsbsim.j3.worst:.2f}} % of peak deviation (criterion 2 %). Most of what is left is about the size of JSBSim's own change when its step is halved; the largest remainder, {{jsbsim.remainder.worst:.2f}} % ({{jsbsim.remainder.where}}), is not attributed to a tested source. The full discrepancy budget is in [`docs/results.md`](docs/results.md#discrepancy-budget). This validates the implementation. It says nothing about how well the tables describe a real F-16.

### Anti-windup rescues a saturating altitude step

![Altitude, pitch command and throttle for a 1,000 ft step, with and without anti-windup](docs/img/c2-windup.png)

A +1,000 ft step saturates the pitch command. With back-calculation anti-windup the altitude overshoots {{control.c2.overshoot:.2f}} % and settles in {{control.c2.settling:.1f}} s; without it, {{control.c2n.overshoot:.2f}} % and {{control.c2n.settling:.1f}} s. The gains come from a linear design (loop margins at least {{control.c5.worst_gm:.0f}} dB and {{control.c5.worst_pm:.0f}}°) and were not tuned on the nonlinear model. Produced by [`f16_ap_c2_alt1000.toml`](scenarios/f16_ap_c2_alt1000.toml) and its [no-anti-windup twin](scenarios/f16_ap_c2_alt1000_nowindup.toml).

The interactive report with every other figure (work–precision, quaternion norm, LN2 and LN3, all JSBSim overlays, C1–C4 and actuator usage) is `out/index.html` after a run.

## Quickstart

Requires CMake ≥ 3.28, a C++20 compiler (GCC ≥ 12, Clang ≥ 15, AppleClang ≥ 15), Ninja, Git and [uv](https://docs.astral.sh/uv/). C++ dependencies (Eigen, GoogleTest, Google Benchmark, toml++) are fetched at configure time with pinned hashes; Python dependencies are locked in `python/uv.lock`. The demo does not need Chrome.

```sh
scripts/demo.sh        # configure, build (release), ctest, uv sync --frozen, run every scenario, report
open out/index.html    # the script prints the path
```

Individual pieces:

```sh
cmake --preset release && cmake --build --preset release -j && ctest --preset release
build/release/flightsim run scenarios/rb_dzhanibekov.toml --out out/rb_dz      # trajectory.csv + metadata.json
build/release/flightsim trim scenarios/f16_trim_502_xcg035.toml --out out/trim.json
build/release/flightsim converge scenarios/studies/core_convergence.toml --out out/convergence
build/release/bench_flightsim
cd python && uv sync --frozen && uv run pytest
cd python && uv run python -m jsbsim_compare.compare --check                  # regenerates the JSBSim references and asserts J1–J3
```

The PNG figures in `docs/img` are exported with `cd python && uv run python -m flightsim_analysis.report --png-dir ../docs/img`; that step needs a Chrome that kaleido can find. `out/index.html` is a single static page with Plotly's JavaScript inlined and works offline in light and dark mode.

## Repository map

```
src/flightsim/
  math/        units, quaternions, frames, Jacobi elliptic functions, interpolation tables
  core/        13-state rigid-body equations of motion, force models, five integrators
  env/         atmospheres (S&L and US 1976), gravity, constant wind and 1-cos gust
  aircraft/    F-16 data loader and the nonlinear aerodynamic, engine and actuator model
  analysis/    trim, central-difference linearization, flight-variable chart, mode identification
  control/     discrete PID with back-calculation anti-windup, altitude and airspeed autopilot, step metrics
  sim/         TOML scenario loaders, runners, CSV and JSON output
  validation/  analytic references (T0–T5, R1–R5) and the convergence harness
apps/          the flightsim command line (run, trim, linearize, converge)
tests/         unit/, validation/ (analytic and published-data checks), alloc/ (zero-allocation test)
bench/         Google Benchmark: integrators and the F-16
scenarios/     every scenario used in a result; control/ holds the autopilot gains; jsbsim/ the matched cases
data/f16/      the S&L tables with citations (f16_sl.toml), reference values, source audit (SOURCES.md)
reference/     generated JSBSim aircraft and planet, and the committed JSBSim result CSVs
python/        flightsim_analysis (run, metrics, plots, report), jsbsim_compare, control_design, tests
docs/          PRD.md (specification), design.md (architecture, decision log), results.md (generated)
scripts/       demo.sh
```

## Validation at a glance

Each row of the validation matrix in [PRD §8.1](docs/PRD.md): what is checked, against which reference, the result, and how strong that reference is.

| Layer | Check | Reference | Result | Strength of the reference |
|---|---|---|---|---|
| Kinematics | T5 rotating body-fixed force, R1 constant spin, R2 symmetric top | Closed forms | Orders in bounds ({{orders.in_bounds}} of {{orders.asserted}} overall); T5 fails for a transposed DCM | Analytic, independent of the code |
| Rigid body | R3 torque-free asymmetric body, H_N and energy conservation | Jacobi elliptic solution (Landau & Lifshitz §37), whose own residuals are checked; conservation laws | Energy drift RK4 {{r3.rk4.dE:.1e}}, RKMK4 {{r3.rkmk4.dE:.1e}}; H_N drift on R4 {{r4.dH:.1e}} | Analytic and conservation laws |
| Integrators | Convergence studies, work–precision | Analytic solutions; for R3 attitude and R4 a fine RK4 reference whose own error is bounded | Orders {{orders.first.min:.2f}}–{{orders.first.max:.2f}} and {{orders.fourth.min:.2f}}–{{orders.fourth.max:.2f}} | Analytic where it exists, otherwise a bounded numerical reference |
| Linearizer | L1 spin-stability eigenvalues | Closed form λ² = −Ω²(I₁−I₂)(I₁−I₃)/(I₂I₃) | Asserted to 1e-6 relative | Analytic |
| Aircraft EOM | B1 work–energy, B2 angular impulse, W1 wind invariance, S1 mirror symmetry | Physical laws | B1 {{checks.b1.residual:.1e}}, B2 {{checks.b2.residual:.1e}}, W1 {{checks.w1.air_relative:.1e}}, S1 {{checks.s1.symmetrized:.1e}} (symmetrized tables) | Physical laws: they involve no published number, and they cannot detect wrong aerodynamic data |
| Aircraft data | Two-source table diff; A1 state-derivative check case | Open reproductions of the S&L model; published numbers as printed by one secondary source | A1 largest relative error {{checks.a1.max_rel:.1e}}; six table disagreements documented in [`resolutions.toml`](data/f16/resolutions.toml) | **Weak to moderate.** The reproductions share one lineage at the root; A1 is a single secondary source (Archimedes) and the book was not available |
| Trim and linear | TR1 published trim, LN1 published linear model, LN2 linear against nonlinear scaling, LN3 step sweep | TR1 and LN1: the same single secondary source; longitudinal matrices there are tutorial-computed. LN2: Taylor's theorem | TR1 {{trim.max_rel_pct:.3f}} %, LN1 eigenvalues {{linear.modes.max_rel_pct:.2f}} %, LN2 ratio {{linear.ln2.ratio:.1f}} | TR1 and LN1: **weak to moderate**, single lineage. LN2 and LN3: analytic |
| Implementation | J1–J3 against JSBSim running a generated copy of the same tables | Independent simulator, same data | J3 worst {{jsbsim.j3.worst:.2f}} % of peak (limit 2 %); J2 trims {{jsbsim.j2.status}} (0.1 %) | Strong for the implementation; says nothing about the data |
| Model form | J4: stock JSBSim `f16` at the same trim | A different model built from TP-1538, Dash-1 and other data | Thrust {{j4.h0.thrust_diff_pct:+.0f}} % and {{j4.h10000.thrust_diff_pct:+.0f}} %, α {{j4.h0.alpha_diff_deg:+.2f}}° and {{j4.h10000.alpha_diff_deg:+.2f}}° (stock minus ours) | Informative only: it sizes the model-form uncertainty |
| Control | C1–C4 step and disturbance metrics, C5 loop margins | The project's own targets (PRD §8.7, tightened at the M7 review) and the design's linear prediction | C1 {{control.c1.overshoot:.2f}} % overshoot; C2 {{control.c2.settling:.1f}} s; C5 margins ≥ {{control.c5.worst_gm:.0f}} dB and {{control.c5.worst_pm:.0f}}° | Self-set targets, checked against a linear model that is itself an approximation |

**Single-lineage caveat.** A1, TR1 and LN1 are compared with the S&L published numbers as reproduced in secondary sources, because the book was not available. [`data/f16/SOURCES.md`](data/f16/SOURCES.md) lists them: the open reproductions descend from one another (the Python and MATLAB AeroBench versions are one lineage; the Archimedes tutorial knows AeroBench), so several citations do not add up to several independent sources, and the longitudinal LN1 matrices are computed by the tutorial rather than printed in the book. Only the Minnesota archive and NASA TP-1538 itself are independent, and they confirm table values, not the check-case outputs. B1, B2, W1, S1, the JSBSim comparison and the analytic cases do not rest on any of this.

## Limitations and future work

Limitations (details in [`docs/design.md`](docs/design.md#5-known-limitations)):

- Flat, non-rotating Earth with NED treated as inertial; constant mass and inertia.
- The F-16 is the reduced S&L model (angle of attack −10° to 45°, subsonic); it is not the full TP-1538 model, and the JSBSim comparison validates the implementation, not the fidelity to a real aircraft.
- A1, TR1 and LN1 rest on secondary sources of a single lineage.
- The JSBSim agreement is limited by JSBSim's own step-halving effect for most signals; the largest remainder ({{jsbsim.remainder.worst:.2f}} %, {{jsbsim.remainder.where}}) is not attributed to a tested source, as the discrepancy budget states.
- The altitude loop has a slow integral tail (about 200 s); the C1 run is 305 s for that reason. The 60°/s elevator rate limit is reached in C1.
- Engine lag is excluded from the JSBSim comparison, and rotor angular momentum is zero on both sides there.
- No turbulence, sensors, noise or state estimation; autopilot is longitudinal only.

Future work: heading hold with coordinated-turn trim, an LQR comparison against the PID cascade, Dryden turbulence, C¹ table interpolation to show RK4's fourth order returning on the aircraft, and access to the book to turn the secondary-source checks into direct transcriptions.

## Documents

- [`docs/PRD.md`](docs/PRD.md): the specification, the validation strategy and every amendment.
- [`docs/design.md`](docs/design.md): architecture, conventions with the test that pins each, derivations, the decision log, known limitations and future work.
- [`docs/results.md`](docs/results.md): every measured number (generated), the benchmark methodology, the J4 appendix and the traceability table.

## License and data attribution

The code is MIT licensed ([`LICENSE`](LICENSE)). The F-16 numbers in `data/f16/f16_sl.toml` come from NASA TP-1538 (Nguyen et al., 1979, a US government report in the public domain), as tabulated by Stevens, Lewis & Johnson, *Aircraft Control and Simulation*, 3rd ed., and as reproduced by the open sources cited in [`data/f16/SOURCES.md`](data/f16/SOURCES.md). Only numbers were transcribed: no book prose or figures and no third-party code were copied, and several of those sources are GPL. JSBSim is a pinned PyPI dependency used from Python and never linked; the stock `f16` of the J4 appendix (GPL, Erik Hofman) is run unmodified and not redistributed.
