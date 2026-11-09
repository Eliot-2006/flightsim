# Design

Companion to `docs/PRD.md` (the specification) and `docs/results.md` (generated numbers). This file explains how the code is organised, why it is built that way, and the mathematics behind the less obvious parts. Everything the PRD specifies is implemented: the rigid-body core, the F-16 model and its consistency laws, trim and linearization, the matched JSBSim comparison and the stock-JSBSim appendix, the actuators and the altitude and airspeed autopilot, benchmarks and the report.

## 1. Architecture

One static library, `flightsim_core`, with folders as modules; a thin CLI; three test executables and a benchmark.

| Module (`src/flightsim/`) | Contents |
|---|---|
| `math/` | units, Hamilton quaternion, frames and Euler angles, Jacobi elliptic functions, interpolation tables |
| `core/` | rigid-body state layout, equations of motion, force models, integrators |
| `env/` | atmosphere models, gravity, constant wind and discrete gust |
| `aircraft/` | F-16 data loader and nonlinear aero/engine/actuator model |
| `validation/` | analytic references (T0–T5, R1–R5), convergence harness |
| `sim/` | TOML scenario loaders (rigid body and aircraft), fixed-step runners including the closed-loop aircraft runner, CSV/JSON output |
| `analysis/` | generic central-difference Jacobian, air-relative flight-variable chart, trim, linearization, mode identification |
| `control/` | discrete PID with back-calculation anti-windup, the altitude/airspeed autopilot, step and disturbance metrics |

Design choices, each of which keeps the dynamics testable in isolation:

- **Pure derivative functions.** `rigid_body_deriv(x, mass_props, force_moment, g_N)` has no state; forces are `std::function` objects that are summed by the caller. The same function is used by the runner, the validation harness and the benchmarks.
- **Integrators are free function templates** over fixed-size `Eigen::Matrix<double, N, 1>`: no heap allocation in a step. `alloc_tests` (`NoAlloc.Rk4StepLoopAllocatesNothing`) replaces global `operator new` and asserts zero allocations over 10⁴ RK4 steps after warmup; `NoAlloc.CounterSeesHeapAllocation` proves the counter is live.
- **Python never computes physics.** `python/flightsim_analysis` runs the CLI (`run`, `trim`, `linearize`, `converge`), reads CSV/JSON, and draws; `python/jsbsim_compare` drives JSBSim and compares; `python/control_design` designs gains from the CLI's linearization. Time-history figures come from scenarios in `scenarios/`.
- **Report pipeline.** `scripts/demo.sh` → build, `ctest`, `uv sync --frozen`, `flightsim_analysis.run` (writes `out/runs/*`, `out/aircraft/*`, `out/control/*`, `out/studies/*`, `out/jsbsim/*` including the stock-f16 trim, `out/checks/*` (the numbers the validation tests print), `out/bench.json` and `out/env.json` (hardware, OS, compiler, flags)) → `flightsim_analysis.report` (writes `out/metrics.json`, `out/index.html`, `docs/results.md`, `README.md`). `docs/results.md` and `README.md` are rendered from `docs/results.template.md` and `docs/README.template.md`; every number is a `{{key}}` from `metrics.json`, and `python/tests/test_report.py` fails if a placeholder is left unfilled.
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

### Decision table

Every entry is dated (the project ran in a single day of working sessions, 2026-11-08), states the decision, the alternatives that were considered and the rationale. The first block is the engineering calls of PRD §1.2; the second the PRD amendments; the third the orchestrator's decisions recorded in `PROGRESS.md`; the fourth the M8 decisions. Detail for the larger items is in the subsections of this file named in the last column.

**Engineering calls (PRD §1.2)**

| Date | Decision | Alternatives | Rationale |
|---|---|---|---|
| 2026-11-08 | Flat, non-rotating Earth; NED treated as inertial | Rotating WGS84 Earth with Coriolis and curvature | Correct for short, low-Mach flights and removes Coriolis and curvature from the core. JSBSim is configured to match; the spike measured 3 to 17 ft of altitude difference over 20 s when it is not (`spikes/jsbsim/FINDINGS.md`). |
| 2026-11-08 | SI units internally; imperial only at data-load and report boundaries | Imperial throughout (S&L, JSBSim) | One unit system in the dynamics; the conversion constants are defined once and tested (`Units.ConstantsAreConsistent`). |
| 2026-11-08 | Translational velocity stored in NED | Body-axis velocity as in S&L and JSBSim | Momentum and energy balances are direct; u, v, w, α, β are derived outputs. |
| 2026-11-08 | Scalar-first Hamilton quaternion `q_NB` (body to NED) | Vector-first or JPL convention; Euler-angle state | Matches Eigen's Hamilton product; the storage-order pitfall has its own test (`Quat.EigenStorageOrderGuard`). |
| 2026-11-08 | Fixed-size Eigen state vectors and integrators as free function templates | Virtual integrator classes; dynamically sized vectors | Zero heap allocation per step (`NoAlloc.Rk4StepLoopAllocatesNothing`) and no class hierarchy. |
| 2026-11-08 | Trim by Levenberg–Marquardt on Eigen with central-difference Jacobians | A general-purpose solver dependency (Ceres, SciPy) | About 50 lines, explainable, no extra dependency (M5 section). |
| 2026-11-08 | Linearize in a 13-state flight-variable chart by central differences | Quaternion coordinates; analytic or complex-step Jacobians | Quaternion coordinates add a spurious zero eigenvalue; analytic derivatives of a table-driven model cost far more than they buy; ε is chosen from the LN3 sweep. |
| 2026-11-08 | TOML scenario files (toml++) | JSON or YAML; command-line flags | Comments document each scenario; Python reads CSV and JSON, never the TOML physics. |
| 2026-11-08 | JSBSim driven from Python, reference CSVs committed | Linking JSBSim into the C++ build | The main build never needs JSBSim; the pinned wheel is a CI job of its own. |
| 2026-11-08 | Controller gains designed on the linear model in Python, verified on the nonlinear model | Tuning on the nonlinear model | Gives margins and a linear-against-nonlinear story; nothing was tuned on the nonlinear model (M7a and M7b sections). |

**PRD amendments**

| Date | Decision | Alternatives | Rationale |
|---|---|---|---|
| 2026-11-08 | §10.3: Euler on R1 has attitude order 2; its first-order behaviour is asserted on the unprojected quaternion-norm error | Keep the PRD's order-1 attitude assertion; drop Euler from R1 | M2 measured order 2: renormalized Euler is an exact rotation about the right axis with a wrong angle (section 3.5). Asserting what is true keeps the test meaningful. |
| 2026-11-08 | §6.3 amendment (M4): for the non-smooth aircraft the 1e-3 reference rule is replaced by RK4 at h_min/64 with a self-convergence estimate, and the fit uses levels whose error is at least 10× that estimate | Keep the rule and claim order 4; skip the aircraft study | With piecewise-linear tables the rule cannot be met; the order is measured and explained, not asserted as 4 (Aircraft budget and convergence amendments). |
| 2026-11-08 | §8.4 amendment (M4): B1 and B2 run 20 s in-domain, not 60 s; a 60 s closed-loop B1 is added at M7 | Widen the table domain; extrapolate | The open-loop doublet excites an unstable mode (+0.13/s) and leaves the α table at 25.46 s; extrapolation would test the tables, not the equations. |
| 2026-11-08 | §8.5 / M5: LN2 runs from a fresh trim at 100 ft, and exactly two inactive-lookup exclusions are allowed, each proved by an in-test assertion | Sea-level run with a loosened ratio bound; no exclusions | At sea level a real lift-loss transient evaluates thrust at h < 0. The [70, 130] bound is unchanged; any other cell change still fails (M5: LN2 active breakpoint guard). |
| 2026-11-08 | §8.6 amendment (spike): JSBSim has no RK4, so its discretization error is bounded by two small steps (1/2400 and 1/4800 s) and its integrators are set explicitly; matched trim needs a dummy zero-power engine; J4 is trim only | Edit a copy of the GPL `f16.xml` to drive its surfaces; run JSBSim at one step | The spike showed the stock FCS overwrites direct surface writes and that FGTrim needs an engine to vary throttle (`spikes/jsbsim/FINDINGS.md`). |
| 2026-11-08 | §8.7 amendment (M7 review): C2 requires anti-windup to be no worse than no anti-windup on overshoot and settling time, in addition to the numeric bounds | Keep the original bounds only | The first implementation met overshoot but never settled; the criterion now encodes what anti-windup is for (M7b: anti-windup tracking times). |

**Orchestrator decisions (`PROGRESS.md`)**

| Date | Decision | Alternatives | Rationale |
|---|---|---|---|
| 2026-11-08 | Commits carry the owner's identity and no co-author or generated-by trailers | Tool-default attribution | The owner's global instruction overrides tool defaults. |
| 2026-11-08 | Python pinned to 3.12 via uv | System Python | One interpreter for every contributor and for CI; the system Pythons on the development Mac are x86 and break `uv` interpreter discovery. |
| 2026-11-08 | CMake minimum raised to 3.28 | Stay at 3.25 and avoid `SYSTEM` in `FetchContent_Declare` | CMake 3.25.1 mis-parses `SYSTEM`; the new minimum was tested in CI and a GCC 13 container. |
| 2026-11-08 | Local sanitizer runs use Homebrew LLVM; CI uses clang-18 | Disable ASan locally | AppleClang 15 ASan aborts at startup on the development Mac. |
| 2026-11-08 | The alloc test is skipped under sanitizers; braces on every `if`/`for` containing a gtest macro | Keep the test in the sanitizer build; accept GCC `-Wdangling-else` | Replacing `operator new` conflicts with the sanitizer allocator; GCC 13 with `-Werror` fails on unbraced gtest macros. A GCC 13 container check runs before every merge. |
| 2026-11-08 | The autopilot operating point is 502 ft/s at 10,000 ft (x_cg 0.35, actuators and engine lag on, 50 Hz) | The PRD's provisional sea-level point | Keeps the gust and step dips away from the ground stop and inside the thrust table (M7a: control operating point). |
| 2026-11-08 | J3 runs at 10,000 ft; J1 and J2 run at sea level and 10,000 ft | J3 from sea level | From sea level the maneuvers leave the thrust table's lower edge, where S&L extrapolates and JSBSim clamps, so the cases would not be comparable. |
| 2026-11-08 | LN1 uses the published xcg = 0.30 eigenvalues; TR1 uses xcg 0.30, 0.35 and 0.38; all labelled single-lineage | Compare LN1 at xcg 0.35 with no reference | The only published linear data found is at 0.30 (M5 section). |
| 2026-11-08 | M4b blocked and resolved: the 60 s open-loop fixture was replaced by a 20 s in-domain run; a 10× reference rule for the aircraft | Continue with the infeasible fixture | See the two amendments above. |
| 2026-11-08 | M5 blocked on the LN2 breakpoint monitor and resolved by the criterion in the LN2 amendment | Drop the monitor | The monitor is what makes the quadratic ratio meaningful. |
| 2026-11-08 | M7 review: back-calculation tracking time tt = T_i for every PI loop, from `design.py`; C2 criterion tightened | Keep tt = 1/ω_c; retune on the nonlinear model | tt = 1/ω_c dragged the altitude integrator negative and the run never settled. The rule is analytic (Åström–Hägglund) and was applied once. |
| 2026-11-08 | M7a blocked twice and resolved: the gains file moved under `scenarios/control/` so the report does not run it as a scenario, and the first design was redone in `design.py` (one cross-feed) | Tune the nonlinear model | The first design missed linear C1/C3 targets; a deterministic search found no feasible plain cascade (M7a: why one cross-feed is needed). |
| 2026-11-08 | The report is one static HTML page with Plotly's JavaScript inlined | A CDN reference | Works offline and cannot drift from a pinned version. |
| 2026-11-08 | Time-history figures use dedicated scenarios in `scenarios/` | Python integration of the plotted trajectories | Everything plotted is produced by the CLI. |
| 2026-11-08 | `docs/results.md` and `README.md` are generated and committed from templates plus a fixture `metrics.json` | Hand-written results | Numbers cannot drift from the code; the fixture keeps the report testable without the C++ build. |
| 2026-11-08 | The integration branch is kept green and PR merges are left to the owner; commits and pushes are the owner's | Orchestrator merges | The owner opened PR #2 and merged PR #1 themselves. |

**M8 decisions**

| Date | Decision | Alternatives | Rationale |
|---|---|---|---|
| 2026-11-08 | J4 reads the stock surface position `fcs/elevator-pos-rad` as δe, and compares thrust in lbf as well as throttle | The stock trim command `fcs/pitch-trim-cmd-norm` as δe; throttle only | The aerodynamic tables read the surface position, and the pitch-trim command passes through a scheduler, PID and kinematic. The throttle scales differ (ours 0.77 is military power, stock 1.0 of a doubled scale is), so thrust is the physical comparison (J4 section). |
| 2026-11-08 | J4 retracts the gear (trimming twice) and lowers the terrain 5,000 ft; the planet stays JSBSim's default | Trim with the file's default gear-down state; build a matched planet for the stock aircraft | The gear adds a large drag term and sea-level ground contact stops FGTrim. J4 is a model-form comparison, so matching the Earth would answer a question the appendix does not ask. |
| 2026-11-08 | Benchmarks report the median CPU time of 7 repetitions of at least 0.1 s each | A single run; the mean | The median is robust to a stray slow repetition on a machine in ordinary use; the coefficient of variation is reported. |
| 2026-11-08 | The F-16 derivative and RK4-step benchmarks start every iteration from the trimmed state | Step on from the previous state | The open-loop aircraft is unstable and leaves the model's domain within seconds of iterations. |
| 2026-11-08 | The 60 s closed-loop benchmark includes its own trim solve | Subtract it; start from a precomputed trim | The trim costs microseconds against a run of milliseconds; including it keeps the benchmark the same code path as the scenario runner. |
| 2026-11-08 | A1, B1, B2, S1, W1 and aircraft-convergence numbers are read from the stdout of the tests that assert them, not duplicated in new code | Add a second computation in C++ or Python; record gtest properties | One computation, one place; a Python test pins the parsed format. |
| 2026-11-08 | The README shows four figures: convergence, Dzhanibekov with H_N, the JSBSim elevator doublet, the C2 anti-windup comparison; `drift.png` is dropped | Keep the drift figure; use C1 | The drift numbers are in the headline table; C2 shows the design insight where C1 only matches the linear prediction. |

## 5. Known limitations

- **Environment.** Flat, non-rotating Earth with NED treated as inertial; constant mass and inertia; uniform gravity unless a model parameter says otherwise. The JSBSim comparison matches a non-rotating spherical planet, and J4 deliberately does not.
- **Reduced S&L data envelope.** The model is the reduced S&L tabulation of TP-1538: α from −10° to 45°, subsonic, narrow β ranges for the control-effectiveness tables, thrust tables over a limited altitude and Mach range. A lookup outside a table extrapolates one cell and is flagged; the validation and comparison runs assert that none occurred, which is why J3 runs at 10,000 ft and B1/B2 run 20 s. Post-stall, spin and landing are outside the model.
- **Single-lineage secondary references.** A1, TR1 and LN1 are compared with the S&L published values as reproduced by secondary sources that descend from one another, because the book was not available; the longitudinal LN1 matrices are tutorial-computed. Several citations do not make several independent sources (`data/f16/SOURCES.md`). B1, B2, W1, S1 and the JSBSim comparison do not depend on them.
- **JSBSim remainder.** For most J3 signals the part of the difference left after the budget is about the size of JSBSim's own change when its step is halved. The largest remainder (the aileron pulse's pitch rate, about 1 %) is not attributed to a tested source: the curvature correction is additive while the aileron pulse swings the attitude by tens of degrees, and small convention differences were not toggled. It is inside the 2 % criterion and stated plainly in `docs/results.md`.
- **Aircraft integration order.** RK4 on the aircraft shows an order near 2, not 4, because the tables are piecewise linear; it is measured and explained, not claimed.
- **Slow altitude-integral tail.** The altitude PI zero at 0.005 rad/s leaves a tail of about 200 s, so C1 runs 305 s and its steady-state criterion is evaluated on the last 10 s.
- **Actuator rate limit.** The 60°/s elevator rate limit is reached at the first controller sample after the C1 command; the linear design has no rate limit, and the comparison shows the effect is brief.
- **Engine lag and rotor excluded from the JSBSim comparison.** Constant throttle and `h_x = 0` on both sides there (JSBSim has no direct equivalent); the production model and the control results include both.
- **No turbulence, sensors or estimation.** Gusts are the deterministic 1-cos profile and a constant wind step; the autopilot sees the true state and is longitudinal only (no heading hold).
- **J4 is trim only.** The stock aircraft's FCS overwrites direct surface writes, so no open-loop response of the stock model is compared, and its planet and mass properties are not matched.
- **Numerics.** Fixed-step integration only; symplectic Euler is exercised on translational oscillators only (§3.4) and is plain Euler on rotation; RKMK4 and Lie–Euler cover attitude only (other states use RK4 and Euler); the zero-allocation test covers the RK4 loop with a quadratic-drag force model; cases without a closed-form attitude (R3 attitude, R4) rely on the fine RK4 reference.
- **Benchmarks.** Numbers in `docs/results.md` are for one machine and one compiler (recorded there, with the methodology) and are reported, not gated.

## 6. Future work

- **Heading hold** (ψ → φ_cmd → aileron with rudder coordination) and coordinated-turn trim (PRD §7.1 "should").
- **LQR** on the same linear model, compared with the PID cascade on the same scenarios and metrics.
- **Dryden / von Kármán turbulence** with a seeded generator, which also needs a sensor and estimation layer to be meaningful.
- **C¹ table interpolation** (for example a monotone cubic or Akima spline) to show RK4's fourth order returning on the aircraft model.
- **Access to the book** (for example through a library): A1, TR1 and LN1 would become direct transcriptions with page citations, and the single-lineage caveat would go. No other change is needed.
- **Closing the JSBSim remainder:** toggling the degree/radian constant and the additive curvature approximation on the aileron case would either explain or bound the remaining 1 %.

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
(The tests now enforce tighter regression bounds; see "Post-M7 review fixes".)
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
`I_k = I* + T (u_sat − u)/tt` when anti-windup is enabled (otherwise I_k = I*). Post-M7: with ki = 0 there is
no integrator (no tracking, and reset stores none), and the pitch-rate loop tracks the measured elevator instead of
`u_sat` (see "Post-M7 review fixes").
Thus integration and back-calculation use forward Euler; the tracking update
requires T < 2 tt for stability. Every generated tt exceeds T by at least 25× (M7b tracking-time rule).
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
recorded for each controller and has no effect with kd = 0. Tracking times follow the
M7b rule (tt = Ti for PI loops; see "M7b: anti-windup tracking times"); the first design used tt = 1/wc. Generated
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
2. *Back-calculation with t_t = 1/ω_c made saturating steps slow* (resolved, see "M7b: anti-windup tracking times" below). With the
   first gains C2 ended 79 ft short after 305 s and C3 did not settle in 80 s, because t_t ≪ T_i.
3. *The elevator rate limit (60°/s) is reached in C1* at the first post-command sample; the linear model ignores it, with negligible
   effect on the C1 comparison.

### M7b: anti-windup tracking times

The first C2 run exposed a design error, not a nonlinear effect. Back-calculation drives the integrator toward its saturated
equilibrium `I* = u_sat − k_p e + t_t k_i e`. With `t_t = 1/ω_c` (altitude 4 s, airspeed 2.5 s, pitch rate 0.125 s) and integral times
`T_i = k_p/k_i` of 200 s, 100 s and 5 s, the term `t_t k_i e` is negligible, so when the proportional term alone exceeded the limit
(C2: θ_cmd; C3: throttle) the integrator was dragged to `u_sat − k_p e`, a large negative value that the slow `k_i` needed hundreds of
seconds to remove: with anti-windup C2 ended 79 ft short after 305 s (no anti-windup: settled in 90.5 s) and C3 never settled in 80 s.
The rule is now `t_t = T_i = 1/zero` for every PI loop (Åström–Hägglund; the PI limit of Åström–Murray's `√(T_i T_d)` with `T_d = 0`), i.e.
altitude 200 s, airspeed 100 s, pitch rate 5 s. With `t_t = T_i` the saturated equilibrium is `I* = u_sat`: while saturated
`dI/dt = (u_sat − I)/T_i`, so the integrator relaxes toward the limit instead of away from it and accumulates less than without
anti-windup, and the loop recovers as soon as P leaves saturation. The P-only pitch loop has no integrator, so its `t_t` keeps the
crossover time. `t_t` does not enter any linear analysis, so C5 margins, linear predictions and all other gains are byte-identical; the
only change to `autopilot_gains.toml` is the three `tt` values. The stability condition `T < 2 t_t` holds with room (every `t_t` is at least 25 T). The rule was chosen from this analysis and applied once; nothing was tuned on the nonlinear model. Result (nonlinear, committed gains):
C2 with anti-windup overshoots 3.08 % (4.61 % without), settles in 21.2 s (90.5 s without), ends 3.46 ft from the command; C3 overshoots 4.22 % and settles in
11.85 s (linear 11.0 s). C1 is unchanged.

### M8: stock JSBSim `f16` (J4)

`python/jsbsim_compare/stock_f16.py` trims the `f16` shipped in the pinned jsbsim wheel (1.3.3) at 502 ft/s, level, at sea level and 10,000 ft, and
`docs/results.md` ("Stock JSBSim F-16 (J4 appendix)") tabulates the result beside our matched-model trim. Four things about the stock file are not
obvious and were found by running it:

- **The gear starts down.** `gear/gear-cmd-norm` defaults to 1 and the file's drag build-up includes a large gear term, so a naive trim describes a
  different aircraft. The script commands the gear up and trims twice: the gear kinematic takes 5 s of FCS time to retract, so the second trim, started
  from the retracted state, is the clean-wing one (`python/tests/test_stock_f16.py` asserts `gear/gear-pos-norm == 0` and a zero gear drag coefficient).
- **Sea-level trim fails on the ground.** The aircraft has bogey and structure contact points, so at 0 ft FGTrim fights ground reactions and reports
  "Trim Failed". `ic/terrain-elevation-ft` is set to −5,000 ft; the atmosphere is a function of altitude above mean sea level and is unchanged.
- **The elevator has three levels.** The pilot input `fcs/elevator-cmd-norm` stays 0. FGTrim varies `fcs/pitch-trim-cmd-norm`, which is summed with the
  stick, clipped, scheduled by angle of attack, fed to a pitch-rate and g-load PID, and then driven through a 0.3 s kinematic to `fcs/elevator-pos-norm`
  and `fcs/elevator-pos-rad` (±0.436 rad). Only the last is the surface position the aerodynamic tables read, so it is the stock equivalent of our δe;
  its sign matches S&L's (the stock pitch-moment table decreases as the position increases).
- **The throttle scale is doubled.** `fcs/throttle-pos-norm = 2 × fcs/throttle-cmd-norm`, and the F100-PW-229 file treats 1.0 as military power and
  1.0 to 2.0 as afterburner. S&L's gearing puts military power at throttle 0.77 (P = 50 %). The numerical throttle values are therefore not comparable;
  thrust in lbf is, and is reported next to the power-lever position as a fraction of military power.

The planet and gravity are JSBSim's defaults (rotating WGS84, gravity model 1), the atmosphere is JSBSim's, and the mass properties are the file's
(20,630 lbf with 3,000 lbf of fuel, CG about 0.018 c̄ forward of the aerodynamic reference point). None of these is matched to ours, by design: the
appendix measures the scale of model-form uncertainty, and no pass/fail applies.

### Post-M7 review fixes

An independent review raised five items; all are fixed without changing any gain, scenario or data value.

1. **P-only loops store no integrator state.** Back-calculation previously ran on every loop, including the P-only pitch loop
   (`ki = 0`). Once that loop saturated (|θ error| > 14.2°), `I += T (u_sat − u)/tt` built an offset that nothing removed: out of
   saturation `u_sat = u` and `ki e = 0`. `reset(u0)` likewise stored `I = u0 − kp e`. Now a `ki = 0` loop never tracks, and reset
   stores no integrator: the first output after reset is still exactly `u0` (bumpless engage), and later outputs follow the plain
   P(D) law. C1–C4 never saturated the pitch loop (C2 peak q_cmd 0.429 against 0.524 rad/s), so this was latent
   (`Pid.ProportionalOnlyLoopKeepsNoBiasAfterSaturationOrReset`).
2. **Pitch-rate anti-windup tracks the actual elevator.** The pitch-rate PI tracked its own ±25° command clamp. But the actuator
   (60°/s, lag 20.2/s) lags a saturated command: in C2 the command reached −25° while the surface reached −14.7°. A new
   `Pid::update(..., applied)` overload back-calculates on `applied − u`. `ApMeasurement::de` carries the applied elevator, which
   is the actuator state when actuators are modelled and the held command for ideal surfaces, sampled with the other
   measurements. In unsaturated flight the actuator lag only shifts the integrator by about −(τ_act + T)/t_t ≈ −1.4 % of each
   command change (τ_act = 1/20.2 s, T = 0.02 s, t_t = 5 s), so the linear design and C5 margins stand and `design.py` and the
   gains are unchanged
   (`Pid.BackCalculationOnAppliedOutputBoundsIntegratorUnderRateLimit`,
   `Autopilot.PitchRateBackCalculationTracksTheMeasuredElevator`). C1–C4 after the change: C1 overshoot 4.36 % (was 4.47 %),
   settling 16.35 s (16.45 s), steady-state 0.335 ft (0.333 ft); C2 3.06 % (3.08 %), 21.25 s (21.2 s), final 3.48 ft (3.46 ft);
   C2 without anti-windup is unchanged (no tracking); C3 4.218 % (4.217 %), 11.85 s; C4 headwind altitude peak 2.365 ft (2.362 ft),
   recovery 136.6 s (136.8 s), gust unchanged to 0.01 ft. Every PRD 8.7 target still passes.
3. **TR1 regression bound.** The PRD acceptance rule (one last printed digit or 0.5 % relative, whichever is larger) is unchanged.
   The test now enforces max(2× observed, one printed digit), which is the printed digit for every TR1 value, since the largest
   observed difference is 0.48 digit. The 0.5 % term had been about 17× the worst observed relative difference (2.8e-4).
4. **LN1 regression bounds** (same rule): eigenvalues max(2.5e-4 relative, digit), where only the short-period imaginary part
   (1.16e-4 relative) exceeds its digit; lateral textbook matrix max(2e-3 relative, digit), worst entry outside its digit 9.23e-4;
   tutorial longitudinal A/B the printed digit alone, because they are rounded to 2–3 significant figures and every entry is within
   0.50 digit.
5. **S&L `adc` pressure.** The routine computes `PS = 1715 ρ T`, a different constant from the 1716.3 in its speed of sound. The
   code had used 1716.3 and now uses 1715. Pressure is output-only (CLI trim JSON, `atmosphere.p`): no force, moment or derivative
   reads it, so no result changes. `Atmosphere.SlAdcMatchesRoutine` reproduces the S&L formulas by construction (offline
   evaluation of the same expressions). It pins the transcription and is not an independent atmosphere check.
