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

### M5: trim and linearization conventions

Trim solves `[Vdot/g, alphadot, qdot]` for alpha, throttle and elevator using
central differences and a three-variable Levenberg–Marquardt system. All three
published level cases converge below the unchanged 1e-10 infinity-norm bound.
TR2 retains its original bounds: the observed 10 s altitude drift at xcg=.35 is
-8.4e-13 ft. At sea level this flags a subpicometer negative-altitude thrust
lookup; the implementation does not clamp the physics or hide that diagnostic.

LN1 uses xcg=.30 exclusively. The full 13-state Jacobian includes power lag;
mode comparisons remove the power state and use `[V,alpha,theta,q]`, holding P
at trim, because the reference has a tabulated engine. The lateral subset is
`[beta,phi,p,r]`. For the tutorial B, throttle uses the algebraic steady-engine
response `B_xu - A_xP B_Pu/A_PP`; elevator columns convert rad to deg, and V
rows/columns convert m/s to ft/s. These operations run in C++, not Python.
The trim result is rejected if its final lookup monitor flags any extrapolation,
including engine Mach/altitude or elevator table limits; CLI trim JSON also
reports trial lookup flags. Climbing-response exports subtract nominal altitude
travel along with horizontal travel, so all exported coordinates are perturbations.
The published longitudinal matrices are tutorial-computed and the lateral
matrix/eigenvalues are textbook values quoted by a secondary source. All share
one Archimedes lineage; multiple citations do not make them independent.

Printed-digit floors from PRD §8.5 matter: tutorial A[alpha,V] prints -0.0003,
where ours is -0.00025372513 (difference +4.62749e-5); the two small B entries
print -0.0021/-0.0022, where ours is -0.0020501361/-0.0021507102. These are
within one last printed digit. Relative bounds remain 1% for textbook values
and 2% for tutorial matrices. No tolerance or source number has been changed.
The small rotational differences remain consistent with the rounded-inertia
finding in A1; there is no calibration of the production inertia.

At xcg=.35 the short-period pair splits into real longitudinal roots, one
unstable. The CLI preserves these as `longitudinal_real_1/2` rather than forcing
a complex-mode classification. The canonical published five-mode comparison
is performed only at xcg=.30.

### M5: LN2 active breakpoint guard (orchestrator ruling)

The original sea-level, rotor-on steps (-0.5 and -0.05 deg, 5 s, RK4 dt=.01)
gave a quadratic ratio of 96.3359912 but flagged four inactive signed-beta
lookups and a brief negative-altitude thrust lookup. The latter is a real
initial lift-loss transient, not an interpolation bug: the large step first
evaluates thrust at h=-4.113491259e-6 m. LN2's altitude is not prescribed by the
plan, so the study now uses a fresh 502 ft/s, xcg=.30 trim at **100 ft** (30.48 m),
inside the same thrust-table cell for the full response. SlAdc, g=data.g_sl,
rotor on, ideal surfaces, both step sizes, the 5 s duration and dt=.01 remain
unchanged. Published TR1/LN1 comparisons still use their specified sea-level
operating points. The CLI exports the separate response operating point and x0.

The orchestrator authorized exactly two exclusions, with assertions in
`F16Linear.LN2QuadraticScaling`:

- DLDA, DLDR, DNDA and DNDR may change cells only while da=dr=0 exactly
  throughout the run. The test asserts ideal surfaces and zero input commands,
  and asserts each run's maximum absolute effective da/dr is exactly zero,
  accumulated at **every RK4 stage**. Each coefficient's output is multiplied
  by zero, so its cell changes cannot act on this response.
- The Cl/Cn sign fold at beta=0 is smooth only inside the same first |beta|
  cell: every loaded alpha-row intercept Cl(alpha,0) and Cn(alpha,0) must be
  exactly zero. The test checks every row of both tables, the unchanged cell
  index and first-beta-cell index. Bilinear interpolation then gives
  `sign(beta)*C(alpha,abs(beta)) = beta*slope(alpha)` across zero. Actual Cl/Cn
  alpha or nonzero-|beta| cell changes are never excluded.

Both step sizes record these events at t=.01 s:

| Table | Cells | Classification and proof |
|---|---|---|
| aero.dlda | 15 → 14 | Inactive: da=dr=0 at every stage |
| aero.dldr | 15 → 14 | Inactive: da=dr=0 at every stage |
| aero.dnda | 15 → 14 | Inactive: da=dr=0 at every stage |
| aero.dndr | 15 → 14 | Inactive: da=dr=0 at every stage |
| aero.cl | 12 → 12 | Smooth beta=0 fold; every alpha-row intercept is zero |
| aero.cn | 12 → 12 | Smooth beta=0 fold; every alpha-row intercept is zero |

The folded Cl/Cn table indices do not actually change; their sign-fold events
are recorded explicitly for audit. The small negative beta is produced by
rotor gyroscopic coupling: `-omega × h_rotor` contributes **q*hx** about body z
when pitch builds after the elevator step, inducing yaw and sideslip.

No other exclusions exist. Any other cx/cz/cm, damping, thrust or Cl/Cn cell
change fails. Regression tests demonstrate that active alpha crossings and a
10,000 ft thrust-altitude breakpoint still fail, nonzero aileron input prevents
the effectiveness exclusion, and one nonzero Cl/Cn intercept at even a remote
alpha row invalidates the fold exclusion. The raw lookup monitor and all
interpolation/domain behavior remain unchanged; classification is confined to
the LN2 diagnostic, and every event includes its table, cells, time and reason
in test output, CLI JSON and the generated report.

At 100 ft the max-error norms are 0.302671106748 and 0.00314170349212, giving
**96.3398065754** within the unchanged [70,130] bound. There are no unexcluded
cell changes and no out-of-domain lookups.

LN3 at xcg=.30 selects epsilon=6e-6. The maximum relative matrix-entry changes
against 6e-5 and 6e-7 are 5.94017e-10 and 1.39313e-10, both below 1e-5. The
1e-2 through 1e-10 sweep exhibits truncation, a flat region and roundoff; each
matrix is normalized by its own maximum absolute entry. The generated report
contains every TR1, LN1 and LN3 difference and the passing LN2 monitor evidence.

### M7a: PID discrete form and engagement

`control/pid.hpp` exposes the Task 7.1 `PidGains` and stateful `Pid` interface.
For sample period T, error e = r − y, and stored integral I, the normal update is
`D_k = (tf D_(k−1) + y_k − y_(k−1))/(tf + T)` (backward Euler derivative filter),
`I* = I_(k−1) + T ki e_k`, `u = kp e_k + I* − kd D_k`,
`u_sat = clamp(u, umin, umax)`, and
`I_k = I* + T (u_sat − u)/tt` when anti-windup is enabled (otherwise I_k = I*).
Thus integration and back-calculation use forward Euler; the tracking update
requires T < 2 tt for stability. Every generated tt exceeds T by at least 6.25×.
Derivative action is on measurement, with its first measurement establishing
history, so a setpoint step cannot kick D. Gains/limits and runtime inputs must
be finite, tf ≥ 0, tt > 0, dt > 0, and umin ≤ umax; invalid inputs throw.

`reset(u0)` defers initialization until the next measurement: clamp u0 to the
output bounds, set D = 0 and I = u0 − kp e, and skip integration for that first
sample. Its output is exactly the requested bounded u0 even with nonzero error.
Construction alone starts with zero integral; explicit reset is the bumpless
engagement operation. The runner should reset altitude to trim theta, pitch to
trim q = 0, pitch-rate to trim elevator, and airspeed to trim throttle. These are
absolute output commands; do not add trim offsets a second time. Disabling
anti-windup only disables tracking, preserving the same saturation and PID law.

### M7a: control operating point and linear model

The orchestrator-selected condition in `scenarios/control/f16_ap_trim.toml` is the S&L
model, V = **502 ft/s**, h = **10,000 ft**, gamma = 0, xcg = 0.35, SlAdc,
g = data.g_sl = 32.17 ft/s², rotor on, actuators **enabled**, engine lag on,
and controller rate **50 Hz**. This replaces the PRD's provisional sea-level-ish
C-scenario point. Trim is alpha = 3.367607324° (theta = alpha), elevator =
−0.656145669°, throttle = 0.168247601, power = 10.925999195 percent.

Reproduce with `cmake --preset release`, `cmake --build --preset release`, then
`cd python && uv run python -m control_design.design`. The script calls
`flightsim linearize` using a temporary copy of the approved scenario with
ideal surfaces, because the current CLI's trim and 13-state linearizer reject
actuator-enabled models. The committed scenario remains actuator-enabled.
This is an analysis-only bypass: Python restores the elevator lag
`20.2/(s+20.2)` from `data/f16/f16_sl.toml`. The six-state longitudinal subset
is `[V,alpha,theta,q,h,P]`, in SI, retaining the **engine power lag** from the
full A/B matrices (no tabulated-engine or frozen-power substitution). The rotor
is enabled during the Jacobian calculation; coupling through lateral states is
neglected by the longitudinal reduction. Task 7.2 must assess that approximation.

Both final command channels (elevator and throttle) include a first-order Padé
approximation to the ZOH half-sample delay, `exp(−sT/2) ≈ (1−sT/4)/(1+sT/4)`.
The cascade is evaluated synchronously within a controller tick, so only the
physical commands incur this hold delay. Elevator position ±25° and rate 60°/s
come from the data file; throttle limits are [0,1]. Position/rate saturation
are inactive in the small-signal model and must be applied in the later runner.
Aileron/rudder commands remain at trim; their dynamics do not enter this
longitudinal design. Altitude output theta_cmd is limited to ±15°; pitch output
q_cmd is limited to ±30°/s as a conservative command bound. The latter is a
controller limit, not an inferred aircraft envelope.

### M7a: loop choices and C5 margins

Pitch rate feeds back q (rad/s) through a **PI** to elevator (rad). The negative
control sign follows S&L's nose-down positive elevator. A zero at 0.2 rad/s
provides low-frequency tracking while retaining phase at the 8 rad/s design
frequency, below the actuator's 20.2 rad/s pole. With all other loops and the
cross-feed closed, crossover is **8.22421 rad/s**, gain-increase margin
**20.8791 dB**, gain-reduction margin **28.0552 dB**, and phase margin
**52.7995°**.

Pitch feeds back theta (rad) through **P** to q_cmd (rad/s). The theta kinematic
integrator already provides steady pitch tracking; another controller integrator
would add avoidable phase lag. Normalize at 2 rad/s after closing the rate loop.
Final crossover is **1.90694 rad/s**, gain-increase margin **18.1372 dB**,
gain-reduction margin **50.0779 dB**, and phase margin **74.4851°**. Its reference
includes the measured-airspeed trim correction described below.

Altitude feeds back h (m, up) through **PI** to theta_cmd (rad). Normalize at
0.25 rad/s after closing the inner loops. A zero at 0.005 rad/s supplies steady
altitude correction without a large integral overshoot; the small slow-pole
residue now lies inside the 2% settling band. Final crossover is
**0.251709 rad/s**, gain-increase margin **20.3648 dB**, and phase margin
**65.0182°**. No finite lower-gain boundary was found. Actual nested-loop
separations are **4.31×** and **7.58×**, both above the required 3× minimum.
The combined pitch reference, including cross-feed, must be limited to ±15°.

Airspeed feeds back air-relative V (m/s) through **PI** to throttle. Normalize at
0.4 rad/s with the altitude cascade closed; the 0.01 rad/s zero gives steady
thrust/drag balance while retaining phase through engine lag. The faster speed
loop limits C1's speed excursion, and the trim-pitch cross-feed prevents that
faster acceleration from lifting the aircraft away from its altitude command.
Final crossover is **0.400000 rad/s**, gain-increase margin **47.3144 dB**,
and phase margin **66.7198°**. No finite lower-gain boundary was found.

### M7a: why one cross-feed is needed

The plain cascade couples a speed change into altitude because the required
level-flight trim alpha (and hence theta) falls as speed increases. Increasing
altitude-loop gain to reject that lift transient worsens altitude-step overshoot;
increasing its integral zero removes the tail but further increases overshoot.
A deterministic differential-evolution exploration of the delayed linear model
(seed 71, 140-member population, up to 250 generations) searched the four
normalization frequencies and three PI zeros. Frequency bounds in rad/s were
q [2,20], pitch [0.2,8], h [0.05,1.5], V [0.03,3]; zero bounds were q [1e-5,20],
h [1e-6,2], V [1e-6,2]. The objective minimized the largest normalized C1/C3
headroom violation with penalties for C5 or separation violations. No feasible
plain-cascade candidate was found: its best compromise still had **9.67% C1
overshoot, 4.83 ft/s speed excursion, and 18.12 ft C3 altitude excursion**,
approximately 1.21 times the requested bounds. This is numerical design evidence,
not a proof of global infeasibility. The much simpler static trim correction
below resolves the coupling without another dynamic controller.

Use exactly one cross-feed, from **measured airspeed** to the pitch reference:
`theta_cmd = clamp(theta_altitude_PI + k_V (V − V_trim), −15°, +15°)`.
Compute k_V from the CLI model's static trim differential: solve
`A delta_x + B[:,(elevator,throttle)] delta_u = 0`, with delta_V = 1 m/s and
delta_h = 0. The h-dot constraint enforces level flight; the power and actuator
states enforce their equilibria. The solution gives
**k_V = d theta_trim/d V = −0.001094945016 rad/(m/s)**, with
**V_trim = 153.0096 m/s**. No fitted compensation or nonlinear trial is used.
Measured speed rather than the stepped speed command avoids an immediate pitch
kick before the aircraft accelerates. No altitude-to-throttle feed-forward is
needed. The cross-feed adds no state, and is included when breaking each loop.
For runner integration, the ±15° bound applies to the **combined** command;
anti-windup tracking must account for that final bound, including the cross-feed.

Synthesis proceeds inner to outer by setting each PI zero and solving the gain
from the return-ratio magnitude at its chosen frequency. Future loops are open
during that loop's synthesis. Final assessment breaks **one error channel at a
time with every other loop closed**, including the airspeed/altitude cross
coupling and measured-speed trim correction. The full longitudinal closed-loop poles must all lie strictly in the
left half plane. Every loop meets C5 (GM ≥ 6 dB, PM ≥ 45°), with an additional
≥6 dB reserve toward any lower-gain stability boundary. Because nested loops
can have both increasing- and decreasing-gain instability, `margins.json`
records both margins and all signed gain boundaries; it does not silently
reinterpret a negative dB boundary as an upper margin. Null lower margins mean
no finite boundary was found on the analysis frequency interval.

Margins use state-space frequency responses over 10⁻⁵…10⁴ rad/s (6,000 logarithmic
samples with interpolated crossing refinement), avoiding spurious crossings
from ill-conditioned high-order polynomial conversion. Crossings beyond
Nyquist are extrapolations of the Padé approximation; the limiting margins
above occur below Nyquist. These are continuous approximations to a sampled
system, not a discrete-time or nonlinear stability certificate. Generated gains
use ten significant digits to remove insignificant compiler roundoff; the
assessment and predictions consume those serialized values. Reruns produce
identical gains, margins and step JSON without timestamps or random choices.

The committed artifact is `scenarios/control/autopilot_gains.toml` (all four seven-field
PidGains tables, rate_hz, and the two crossfeed constants). Derivative gains are zero; tf = T = 0.02 s is
recorded for each controller and has no effect with kd = 0. Tracking tt = 1/wc
sets recovery on each loop's nominal response time. Generated
`build/control_design/margins.json` contains the margin report and
`build/control_design/linear_steps.json` contains C1 (+100 ft = 30.48 m) and
C3 (+50 ft/s = 15.24 m/s) responses on a 0.1 s grid through 300 s. These are
**perturbations** from trim, with steps at t = 0: add the trim h and V for
absolute overlays, and shift time to the runner's command event. No nonlinear
run is used in design or acceptance tests.

The final delayed linear predictions meet the orchestrator's headroom targets:

| Case | Metric | Prediction | Linear design target |
|---|---|---:|---:|
| C1 +100 ft | Altitude overshoot | 4.3021% | ≤8% |
| C1 | 2% settling | 16.3 s | ≤20 s |
| C1 | Peak absolute speed excursion | 2.9582 ft/s | ≤4 ft/s |
| C3 +50 ft/s | Speed overshoot | 4.2835% | ≤8% |
| C3 | Peak absolute altitude excursion | 4.1807 ft | ≤15 ft |

C3 also settles within 2% in 11.0 s. The previous 50.4 s C1 tail is removed
from the settling band; no known linear target miss is deferred to Task 7.2.
`generate()` refuses to emit a design that fails these targets, C5, pole
stability, or bandwidth separation. Tests independently check the predicted
traces against the specified bounds as well as reproducibility. Nonlinear
verification, rate/position saturation and exact sampled-system behavior remain
Task 7.2 checks; these linear results do not claim those checks have passed.

### M7b: closed-loop runner, autopilot and nonlinear verification

**Trim with actuators.** `trim_level` solves on an ideal-surface copy of the model, then returns the 17-state vector
with the elevator/aileron/rudder states equal to the trimmed commands (`x17_from_x13` sets them from `u`) and the
power state at its commanded steady value, so the actuator and power derivatives are zero to roundoff
(`Trim.ActuatorEnabledTrimIsEquilibriumForActuatorAndPowerStates`). `flightsim trim` therefore accepts actuator-enabled
scenarios; `linearize` still requires ideal surfaces, and `design.py` keeps its documented bypass.

**Autopilot** (`control/autopilot.*`) is the M7a structure and gains verbatim: altitude PI → θ_cmd, plus the single cross-feed
`k_V (V − V_ref)` on *measured air-relative* V, combined limit ±15° (the altitude gains' own bounds); pitch P → q_cmd; pitch-rate
PI → elevator; airspeed PI → throttle. Measurements are x13 values (h, V, θ, q). The altitude PI is called with per-sample limits
`[lo − xf, hi − xf]` (clamped into its own bounds), so its back-calculation tracks the *combined* bound, as M7a requires; this is
the only change to `Pid` (an `update` overload taking limits; the 3-argument form forwards to it). `engage` resets the PIDs to
trim: altitude → θ_trim − k_V(V − V_ref), pitch → 0, pitch rate → trim elevator, airspeed → trim throttle, so the first sample
outputs the trim command exactly (`EngageIsBumplessAndHoldsTrimAtZeroError`). Unit tests pin the signs and the cross-feed with
hand-computed values from the committed gains.

**Runner.** `run_aircraft` calls the controller before the integration step that starts on a sample instant
(`k % steps_per_sample == 0`) and holds `u` until the next sample. `dt_s` must divide `1/rate_hz` (load error), and
`autopilot.rate_hz` must equal the rate the gains were designed for (load error). Commands are absolute setpoints
(`[[commands]] kind, t_s, value_ft|value_m|value_fps|value_mps`), picked up at the first controller sample at or after `t_s`; before
a command the setpoint is the trimmed value. `[[inputs]]` and `[autopilot]` are mutually exclusive. The headwind change is a
time-aligned constant-wind step, `[[wind.steps]] t_s, delta_mps`: the wind in effect over a step is the base constant plus all steps
with `t_s ≤ t` (`constant_wind_at`), and `t_s` is a multiple of `dt_s`, so no step straddles it.

**Metrics** (`control/metrics.*`, pure functions; `control_report` assembles them): overshoot = max excursion past the final value in
the step direction, % of |step|; 2 % settling = time from the command to the first sample after which |y − y_final| stays within
2 % of |step| (null if never); steady-state error = max |y − y_final| over the last 10 s of the record; the other loop's peak error;
actuator peaks over every integration step. Disturbance metrics use windows from one disturbance to the next: peak |y − command| and
the time to stay within 2 % of that peak. Metrics are written to `metadata.json` under `control`; a time that was not reached is
omitted and flagged `*_reached = false` (JSON has no NaN).

**B1 closed loop.** `F16AutopilotB1.ClosedLoopWorkEnergyBalance60s` integrates `[x17 | W_nc]` through 60 s of C1 with the controller
in the loop. The balance holds per step because the controller is outside the derivative: inputs are constants within a step
(sample and hold, `dt` divides the period), so each step integrates the same smooth ODE as the open-loop run, and
d(E − W_nc)/dt = 0 holds for that ODE under any input; the controller only changes the constants between steps, while E and W_nc
are functions of the continuous state. The test's hand-rolled loop is checked against `run_aircraft` (same final state).

**Nonlinear results and what they say about the design.** C1 matches the linear prediction closely (see results.md). Three items
for the M7 tuning review, all diagnosed, none changed (gains go through `design.py`):

1. *C1 steady state is horizon-limited.* The altitude PI zero at 0.005 rad/s leaves a tail of about 200 s; the linear model itself has
   |e| = 1.09 ft 55 s after the command and 0.32 ft at 300 s, and the nonlinear run agrees to 0.001 m. The C1 scenario runs the 300 s
   horizon of `linear_steps.json` (chosen from the linear prediction) and the 1 ft target is evaluated on its last 10 s.
2. *Back-calculation with t_t = 1/ω_c makes saturating steps slow.* When P alone exceeds an actuator limit (C2: θ_cmd; C3: throttle,
   k_p·ΔV = 1.09 > range), back-calculation drives the integrator to u_sat − k_p e, a large negative value that the loop's slow
   k_i (altitude 8.6e-6, airspeed 7.2e-4) takes hundreds of seconds to undo. C2 with anti-windup overshoots less (asserted) but ends
   far short of the command after 305 s; without it the response settles in about 90 s. C3's airspeed never settles within 2 % in
   80 s, whereas a scratch run (not committed) with anti-windup off reproduced the linear prediction (overshoot ≈ 4.6 %, settling
   ≈ 12.7 s). The tracking time t_t is the lever; PRD targets C2 (overshoot ≤ 15 %) and C3 (overshoot ≤ 10 %, altitude ≤ 20 ft) are met.
3. *The elevator rate limit (60°/s) is reached in C1* at the first post-command sample; the linear model ignores it, with negligible
   effect on the C1 comparison.
