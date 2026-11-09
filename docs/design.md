# Design

Companion to `docs/PRD.md` (the specification) and `docs/results.md` (generated numbers). This file explains how the code is organised, why it is built that way, and the mathematics behind the less obvious parts. The F-16 model and validation are implemented; trim, aircraft linearization, JSBSim comparison and autopilot follow in M5–M7.

## 1. Architecture

One static library, `flightsim_core`, with folders as modules; a thin CLI; three test executables and a benchmark.

| Module (`src/flightsim/`) | Contents |
|---|---|
| `math/` | units, Hamilton quaternion, frames and Euler angles, Jacobi elliptic functions, interpolation tables |
| `core/` | rigid-body state layout, equations of motion, force models, integrators |
| `env/` | atmosphere models, gravity, constant wind and discrete gust |
| `aircraft/` | F-16 data loader and nonlinear aero/engine/actuator model |
| `validation/` | analytic references (T0–T5, R1–R5), convergence harness |
| `sim/` | TOML scenario loader, runner, CSV/JSON output |
| `analysis/` | generic central-difference Jacobian, air-relative flight-variable chart |

Design choices, each of which keeps the dynamics testable in isolation:

- **Pure derivative functions.** `rigid_body_deriv(x, mass_props, force_moment, g_N)` has no state; forces are `std::function` objects that are summed by the caller. The same function is used by the runner, the validation harness and the benchmarks.
- **Integrators are free function templates** over fixed-size `Eigen::Matrix<double, N, 1>`: no heap allocation in a step. `alloc_tests` (`NoAlloc.Rk4StepLoopAllocatesNothing`) replaces global `operator new` and asserts zero allocations over 10⁴ RK4 steps after warmup; `NoAlloc.CounterSeesHeapAllocation` proves the counter is live.
- **Python never computes physics.** `python/flightsim_analysis` runs the CLI (`run`, `converge`), reads CSV/JSON, and draws. Time-history figures come from scenarios in `scenarios/`.
- **Report pipeline.** `scripts/demo.sh` → build, `ctest`, `uv sync --frozen`, `flightsim_analysis.run` (writes `out/runs/*`, `out/studies/*`, `out/bench.json`) → `flightsim_analysis.report` (writes `out/metrics.json`, `out/index.html`, `docs/results.md`). `docs/results.md` is rendered from `docs/results.template.md`; every number is a `{{key}}` from `metrics.json`, and `python/tests/test_report.py` fails if a placeholder is left unfilled.
- **Extending the report.** A new figure is one builder plus one `Figure(...)` entry in `plots.py`; a new result is keys in `metrics.py` plus text in the template. The page is a single static file with Plotly JS inlined (no CDN, no server).
- **Determinism.** `-ffp-contract=off` on all targets; `Runner.DeterministicBitwise` asserts two runs write identical bytes.

## 2. Conventions (PRD §4) and the test that pins each

| Convention | Pinned by |
|---|---|
| Frames NED / FRD; `C_NB` maps body to NED | `Quat.RotateBasis90` (ψ=90°: body x → East; θ=90°: body x → Up; φ=90°: body y → Down), `RigidBody.BodyForceRotatedToNed`, `Translational.T5DirectionOfRotation` (a transposed DCM diverges by > 0.1) |
| Euler 3-2-1, `C_NB = R_z R_y R_x` | `Quat.EulerSequence321MatchesElementaryRotations`, `Quat.EulerRoundTrip` |
| Gimbal lock at θ = ±90° defined (φ = 0) | `Quat.GimbalLockDefined` |
| Euler-rate matrix | `Frames.EulerRateMatrixMatchesFiniteDifference` |
| α = atan2(w, u), β = asin(v/V), V → 0 guarded | `Frames.AlphaPositiveForPositiveW`, `Frames.BetaPositiveForPositiveV`, `Frames.ZeroAirspeedDefined` |
| Hamilton, scalar-first, `q_NB` | `Quat.HamiltonProductMatchesHandComputed`, `Quat.ConjugateInverse`, `Quat.RotateMatchesSandwichProduct` |
| Kinematics `q̇ = ½ q ⊗ [0, ω_B]` | `RigidBody.QuatKinematicsMatchesFormula`; `Rotational.R2Rk4Order4` |
| Eigen `coeffs()` is `(x, y, z, w)` | `Quat.EigenStorageOrderGuard` |
| q and −q equal; error angle sign-invariant | `Quat.ErrorAngleSignInvariant` |
| Renormalize after a full Euler/RK4 step; Lie methods do not | `Integrators.RenormalizeLeavesOthers`, `Rotational.NormDriftWithoutRenormalization`, `LieIntegrators.LieEulerPreservesNorm` |
| SI units, exact imperial constants | `Units.ConstantsAreConsistent` |
| 13-state layout `[p_N, v_N, q_NB, ω_B]` | `RigidBody.StateLayoutAndAccessors` |
| Gravity as a model parameter | `RigidBody.GravityOnlyAcceleration` |
| Scenario `dt` divides run, output interval | `Scenario.RejectsDtNotDividingTFinal` |
| Non-finite state stops the run | `Runner.NonFiniteStops` |

The aircraft layout and wind-relative chart are pinned by `FlightVars.RoundTripX13`
and `F16Consistency.W1GalileanInvariance`; the S&L control signs are pinned by
`F16Model.PositiveElevatorNoseDown`, `PositiveAileronLeftRoll` and
`PositiveRudderLeftYaw`.

## 3. Derivations

### 3.1 Equations of motion

With state `x = [p_N, v_N, q_NB, ω_B]`, mass `m`, inertia `I` about the CG, rotor momentum `h`, non-gravitational body force `F_B` and moment `M_B`:

```
ṗ = v
v̇ = C_NB F_B / m + g_N          (F_N added directly when a force is specified in NED)
q̇ = ½ q ⊗ [0, ω_B]
ω̇ = I⁻¹ ( M_B − ω × (I ω + h) )
```

The first is kinematics, the second Newton's law in the (inertial) NED frame, the third quaternion kinematics, the fourth Euler's rotational equation in the body frame. Torque-free, zero-rotor motion conserves `E = ½ ωᵀIω` and `H_N = C_NB I ω`, which is what the R-series tests check. `RigidBody.EulerEquationsGyroscopic` and `RigidBody.RotorTermIncluded` fix the signs of the gyroscopic and rotor terms.

### 3.2 Quaternion kinematics

For a body rotating with angular velocity ω_B (body components), an infinitesimal rotation over `dt` composes on the **right**: `q(t+dt) = q(t) ⊗ [1, ½ ω_B dt]`. Hence `q̇ = ½ q ⊗ [0, ω_B]`. Writing it `½ [0, ω] ⊗ q` instead would mean ω is expressed in N and would silently rotate the wrong way for body rates; R2 (non-trivial attitude, tilted ω) and `Translational.T5DirectionOfRotation` exercise the body-rate convention with a generic attitude; a single-axis spin about the initial axis (R1) cannot distinguish the two.

### 3.3 RKMK4 and the +½ sign

Write `q(t) = q_k ⊗ exp(½ Θ(t))`, Θ ∈ ℝ³, Θ(0) = 0 (body coordinates). Differentiating gives `Θ̇ = dexp⁻¹_{−Θ}(ω_B(t))`: the update multiplies on the right, so the inverse right Jacobian of the exponential map is evaluated at −Θ. Its series truncated at the term needed for order 4 is

```
Θ̇ = ω + ½ Θ × ω + (1/12) Θ × (Θ × ω)
```

The `+½` follows from the argument being **−Θ**: `dexp⁻¹_{Θ}(ω) = ω − ½ Θ×ω + (1/12) Θ×(Θ×ω) + …`, and replacing Θ by −Θ flips the odd term. The code (`rkmk4_step` in `core/integrators.hpp`) applies classical RK4 to Θ and to the vector-space states, then sets `q_{k+1} = q_k ⊗ exp(½ Θ_{k+1})`. With `−½` the method degrades to order 2 on R2/R3 while R1 (Θ ∥ ω, so the cross term vanishes) still passes; `LieIntegrators.Rkmk4SignRegression` compares one step against a tiny-step RK4 reference and requires the `−½` variant, computed in the test, to be at least 100× worse. The norm is preserved to roundoff with no projection.

### 3.4 Symplectic Euler: scope

`v⁺ = v + h a(p)`, `p⁺ = p + h v⁺` is symplectic for separable Hamiltonians `½|v|²/m + U(p)` (T4, a harmonic oscillator): energy error stays bounded and does not trend (`LieIntegrators.SymplecticEulerOscillatorEnergyBounded`; Euler's grows, `Translational.T4EulerEnergyGrows`). It is **not** applied as a symplectic method to rotation (the attitude/ω flow is Lie–Poisson, not canonical in these coordinates) or to velocity-dependent drag (T2/T3): in those cases the code uses plain Euler for the other states, and the PRD §10.3 table asserts symplectic Euler on T4 only. Its observed order on T4 is still 1.

### 3.5 Euler on R1 has attitude order 2 (PRD §10.3 amendment)

R1 is constant ω about a principal axis. For Euler with renormalization each step is `q⁺ = q + (h/2) q ⊗ [0, ω] = q ⊗ [1, (h/2) ω]`. After normalization this is an **exact rotation about ω̂** (the axis is right), through the angle `2·atan(h|ω|/2) = h|ω| − (h|ω|)³/12 + …`. Only the angle is wrong, by O(h³) per step; over T/h steps the error is O(h²). So attitude error has order 2, not 1. The first-order behaviour of Euler on R1 lives in the **unprojected** quaternion norm: each step multiplies the norm by `√(1 + (h|ω|/2)²)`, so after T/h steps the norm error is ≈ T|ω|² h / 8, order 1. The study therefore asserts Euler on R1 with `renormalize = false` on `err_qnorm` (`Rotational.R1EulerOrder1`) and reports the attitude order without asserting it. RK4 on R1 keeps the attitude assertion.

### 3.6 Convergence methodology

Implemented in `validation/convergence.*`; normative text in PRD §6.3.

- Fixed final time, `h_k = h₀/2ᵏ`, eight levels.
- Reference: the analytic solution where one exists. Otherwise an RK4 run at `h_K/16` regardless of the method under test (a same-method reference biases first-order orders upward), refined until its own error, estimated by step halving, is below 10⁻³ of the finest-level error of the method under test.
- Metrics: position, velocity and rate Euclidean norms; attitude via the sign-invariant error angle; relative energy and `H_N` drift.
- Observed order: least-squares slope of log error against log h over the asymptotic window, i.e. levels whose error exceeds 10³·ε·scale and whose successive local orders agree within 10 %. The asserted bounds are those of PRD §10.3 ([0.95, 1.05] and [3.9, 4.1]); `docs/results.md` lists every observed order.
- Cases where a method is exact or roundoff-limited (RKMK4 on R1) assert the absolute error instead (`LieIntegrators.Rkmk4ExactForConstantOmega`).
- Output: `convergence.csv`, `orders.json`, `norm_conservation.csv`; the report plots error against h with reference slopes and against derivative evaluations (work–precision).

## 4. Decision log

### M4 wind profile

`WindModel` adds a constant NED wind and an optional MIL-F-8785C 1-cos discrete
gust, applied as a full rise and fall over its duration. The peak is the specified
amplitude along the normalized NED direction at the midpoint; endpoints are zero.
This is a frozen time profile: spatial gradients and rotational gusts are neglected
(PRD §5.4). Event alignment belongs to the scenario caller.

### M4 F-16 dynamics and chart

`F16Model` assembles the S&L forces and moments using the SI-converted loader
tables. Static Cl/Cn use odd folding in beta; effectiveness tables use signed beta.
CG transfer uses total Cy/Cz after damping. The control divisors are 20°/30°,
independent of the actuator position limits. Mass is the reciprocal of S&L RM;
the FRD inertia tensor has off-diagonal -Ixz and engine momentum along positive
body x. The generic rigid-body EOM owns the gyroscopic term.

`flight_vars` maps between NED ground velocity and the air-relative S&L chart.
For steady wind, body air acceleration is `C_BN * vdot_N - omega.cross(v_air_B)`;
the chart derivatives use that acceleration, including in nonzero constant wind.
Below 1e-9 m/s the existing air-data convention defines alpha/beta as zero and
the model suppresses nondimensional damping rates. At rest, aerodynamic force
vanishes and the NED state can accelerate under thrust/gravity. The singular
V/alpha/beta chart derivatives are defined as zero there; use the 17-state model
to depart from rest. Euler pitch and sideslip chart poles throw `domain_error`.

Throttle is clamped to its data-file limits. With actuators enabled, actual
surface states determine the force and lag/rate-limited derivatives drive them
toward position-limited commands. With actuators bypassed, commands determine
force directly and surface-state derivatives are zero. `Eval.surfaces` reports
effective throttle followed by elevator, aileron and rudder (radians).

A1 is a secondary, single-lineage check (Archimedes). The full inertia gives
small differences from the published rounded c1…c9 rotational equations; the
test reproduces those rounded equations diagnostically without changing the
production EOM. The investigation is recorded under Validation findings below.

From PRD §1.2 (engineering calls made with the specification):

| Call | Rationale |
|---|---|
| Flat, non-rotating Earth; NED treated as inertial | Correct for short-duration, low-Mach flight; removes Coriolis/curvature. JSBSim is configured to match. |
| SI units internally; imperial only at data-load/report boundaries | One unit system in the dynamics; conversion is tested. |
| Translational velocity stored in NED | Direct momentum/energy balances; u, v, w, α, β are derived outputs. |
| Scalar-first Hamilton quaternion `q_NB` (body → NED) | Matches Eigen's Hamilton product; documented once and tested. |
| Fixed-size Eigen state vectors; integrators as free function templates | Zero heap allocation per step; no class hierarchy. |
| Trim by damped Newton / Levenberg–Marquardt on Eigen | Small, explainable, no extra dependency. *(M5)* |
| Linearization by central differences in a 13-state flight-variable chart | Avoids the redundant-quaternion zero eigenvalue. *(M5; generic `jacobian_central` exists)* |
| TOML scenario files (toml++) | Human-readable, comments document scenarios; Python reads CSV/JSON. |
| JSBSim driven from Python, reference CSVs committed | No C++ dependency on JSBSim. *(M6)* |
| Controller gains designed on the linear model in Python | Classic workflow, margins plus nonlinear verification. *(M7)* |

From `PROGRESS.md` and M2/M3 work:

| Date | Decision |
|---|---|
| 2026-11-08 | Commits carry the owner's identity and no co-author or generated-by trailers. |
| 2026-11-08 | Python pinned to 3.12 via uv. |
| 2026-11-08 | LN1 uses the xcg = 0.30 published eigenvalues; TR1 uses xcg 0.30/0.35/0.38 level trims; all labelled single-lineage per PRD §8.5. *(M5)* |
| 2026-11-08 | PRD §10.3 amended: Euler on R1 has attitude order 2; asserted on norm error instead (§3.5). |
| 2026-11-08 | PRD §8.6 amended per the JSBSim spike: no-RK4 JSBSim reference via two small steps; J4 trim-only. *(M6)* |
| 2026-11-08 | Report is one static HTML page with Plotly JS inlined (≈ 5 MB), not a CDN reference: works offline and cannot drift from a pinned version. |
| 2026-11-08 | Time-history figures use dedicated scenarios in `scenarios/` (R3 at h = 0.02 s per integrator, T4 at h = 0.05 s) so everything is produced by the CLI. |
| 2026-11-08 | `docs/results.md` is generated and committed; the template and a fixture `metrics.json` keep it testable without the C++ build. |
| 2026-11-08 | `bench_smoke.cpp` removed; `bench_integrators.cpp` supersedes it. |

## 5. Known limitations

- Aircraft trim, linearization, JSBSim comparison and autopilot are *coming in M5–M7*.
- Flat non-rotating Earth, constant mass and inertia, uniform gravity unless a model parameter says otherwise.
- Fixed-step integration only; no adaptive stepping or dense output.
- Symplectic Euler is exercised on translational oscillators only (§3.4); on rotation it is plain Euler.
- RKMK4 and Lie–Euler cover attitude only; other states use RK4 and Euler respectively.
- Smooth-case observed orders use eight step levels with reference-accuracy checks; cases without a closed-form attitude (R3 attitude, R4) rely on the fine RK4 reference. The non-smooth aircraft study uses the amended five-level procedure below.
- Zero-allocation is asserted for the RK4 loop with a quadratic-drag force model; other integrators and force combinations are not separately tested.
- Benchmark numbers in `docs/results.md` are for one machine and one compiler (recorded there); they are reported, not gated.

## Validation findings

### A1: rounded rotational inertia coefficients

A1 uses the S&L Table 3.5-2 values reproduced by Archimedes, a **single
reference lineage** (see `data/f16/reference_values.toml` for citations).
`F16Model.A1CheckCase` computes in SI, converts back to feet/percent for
comparison, and prints every component's relative error. All components retain
the binding 1e-3 relative tolerance.

The two errors above the 1e-4 investigation threshold are pdot
(1.99905157675e-4) and rdot (1.03599879618e-4). Production derives rotational
coupling from the full inertia tensor, with off-diagonal -Ixz. The printed
reference instead uses the rounded S&L c1…c9 combinations. Evaluating those
rounded coefficients diagnostically yields
`[12.6267851355, 0.964966917861, 0.580975825351] rad/s²`, matching the printed
`[12.62679, 0.9649671, 0.5809759]` to within 4e-7 relative. The test checks this
explanation and exempts only pdot/rdot from the investigation failure, retaining
the 1e-3 bound for both. It does not alter production inertia or the generic EOM.
The calculation uses rtod=57.29578, the distinct 57.3 Cz constant, mass=1/RM,
g=32.17 ft/s², SlAdc and the negative-beta odd fold.

### Aircraft budget and convergence amendments

The prescribed open-loop doublets excite an unstable longitudinal mode
(+0.13257/s near the initial point), leaving the alpha table domain at 25.46 s.
The orchestrator amended B1, B2 and the aircraft convergence run to **20 s**,
with the doublets and domain assertions unchanged. B1/B2 retain their 1e-6
relative residual bound at dt=0.01. A 60 s closed-loop budget check belongs to M7.

For the non-smooth aircraft convergence study, the amended reference is RK4 at
`h_min/64`, with error estimated from its difference against `h_min/32`. Each
metric's fit excludes levels with error below 10 times that reference-error
estimate; the test prints the estimate and the actual fitted levels. It fits
the least-squares log-error/log-step slope over all reference-resolved levels
and retains the [0.9, 4.2] order bounds. Piecewise-linear tables reduce smoothness
at cell crossings, so aircraft order is measured rather than assumed to be four.
The smooth-case harness, its asymptotic-window selection and its 1e-3
reference-accuracy rule are unchanged.
