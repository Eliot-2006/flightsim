# PRD — 6-DOF Flight Dynamics Simulator

| | |
|---|---|
| Status | Draft for approval |
| Date | 2026-11-08 |
| Owner | Eliot Boda |
| Language / build | C++20, CMake ≥ 3.25 |
| Next step after approval | Implementation plan per milestone (superpowers:writing-plans) |

---

## 1. Intent

Build a six-degree-of-freedom rigid-body flight dynamics simulator whose value is **demonstrated correctness**, not visuals. Every physical claim in the README must be backed by a test or a reproducible measurement in this repository. The project should hold up under questioning in aerospace, robotics, and simulation-engineering interviews: every convention, approximation, and discrepancy is written down with a reason.

The project has two layers, kept strictly separate:

1. **Generic rigid-body engine.** 6-DOF equations of motion, quaternion attitude, pluggable forces/moments, multiple integrators. Validated against closed-form solutions and conservation laws.
2. **Aircraft flight dynamics model.** F-16 nonlinear model from published NASA wind-tunnel data (as tabulated by Stevens, Lewis & Johnson), atmosphere, engine, actuators, trim, linearization, and an altitude/airspeed autopilot. Validated against published check cases, published trim/linear results, internal consistency laws, and a matched-model JSBSim comparison.

### 1.1 Decisions made with the owner

| Decision | Choice |
|---|---|
| Aircraft model | F-16, Stevens, Lewis & Johnson (*Aircraft Control and Simulation*, 3rd ed., Wiley 2016, "S&L") Appendix A model; data originates from NASA TP-1538 (Nguyen et al., 1979, public domain) |
| Required control scope | Altitude hold + airspeed hold (cascaded PID, anti-windup, actuator limits). Heading hold is stretch. |
| Depth extras in scope | Lie-group (exp-map, RKMK4) attitude integrator; symplectic Euler; appendix comparing against JSBSim's stock F-16 |
| Schedule | No fixed deadline. Milestones are gated by acceptance criteria; effort estimates are given for planning only. |
| Book access | Owner does **not** have S&L. Tables and published reference values come from open reproductions of the S&L model, cross-checked between ≥ 2 independent sources, with NASA TP-1538 as primary citation (§5.3, §8.5). |
| Data in public repo | Commit the numeric tables in `data/f16/f16_sl.toml`, citing NASA TP-1538 as primary source and S&L as the tabulation. No book prose, figures, or third-party code copied. |
| Review gates | Owner signs off every "manual review" gate. |

### 1.2 Engineering calls made in this PRD (owner may override)

| Call | Rationale |
|---|---|
| Flat, non-rotating Earth; NED treated as inertial | Correct for short-duration, low-Mach flight; removes Coriolis/curvature from the core model. JSBSim is configured to match (§9.3). |
| SI units internally; imperial only at data-load / reporting boundaries | One unit system in the dynamics. S&L and JSBSim data are imperial; conversion is tested. |
| Translational velocity stored in NED, not body axes | Makes linear-momentum/energy balances direct; aircraft quantities (u,v,w,α,β) are derived outputs. |
| Scalar-first Hamilton quaternion `q_NB` (body→NED) | Industry-standard; matches Eigen's Hamilton product; documented once and tested. |
| Fixed-size Eigen state vectors, integrators as free function templates | Zero heap allocation per step, no class hierarchy, easy to reason about. |
| Trim via damped Newton / Levenberg–Marquardt written on Eigen | ~50 lines, explainable, no extra dependency. |
| Linearization by central finite differences in a 13-state flight-variable chart | Avoids the spurious zero eigenvalue of the redundant quaternion; matches S&L state definitions. |
| TOML scenario files (toml++) | Human-readable, supports comments to document scenarios; Python reads with stdlib `tomllib`. |
| JSBSim driven from Python (`pip install jsbsim`, pinned), reference CSVs committed | No C++ link dependency on JSBSim; main build never needs it. |
| Controller gains designed on the linear model in Python (`python-control`), verified on the nonlinear model | Classic workflow; produces margins and a "linear vs nonlinear" story. |

---

## 2. Critical review of the original brief

Findings that shaped this PRD:

1. **Some proposed validations are tautologies.** Linear momentum conservation with zero force on a single body is exactly preserved by every integrator; constant acceleration is integrated exactly by RK4. Both are kept only as smoke tests and are *not* cited as evidence of correctness.
2. **The strong rotational tests are nonlinear.** Torque-free motion of an asymmetric body (Jacobi elliptic solution), the axisymmetric top (closed-form attitude), and conservation of the **inertial-frame** angular-momentum vector (which tests ω *and* attitude integration together).
3. **"Energy conservation" does not apply to an aircraft.** Aircraft get **work–energy and angular-impulse balances** instead: ΔE (incl. rotational KE) must equal ∫(non-gravitational power)dt; ΔH_N must equal ∫M_N dt. These validate the EOM and attitude kinematics independently of aerodynamic data correctness (transform *direction* is covered separately by a closed-form test).
4. **RK4 will not show 4th order on the aircraft model.** Linearly interpolated tables, saturations, and input switching are non-smooth. The PRD requires measuring and explaining the reduced order, not claiming 4.
5. **A JSBSim comparison is only meaningful if the models match.** JSBSim's stock F-16 appears to use a different (larger) TP-1538 data set and includes a flight-control system. The required comparison uses a JSBSim aircraft file **generated from the same data file** our C++ code loads, with Earth rotation, gravity model, inertia sign convention, table extrapolation, and atmosphere aligned. That comparison validates the *implementation*; aero *fidelity* rests on data provenance.
6. **The S&L model has its own atmosphere and gravity constant.** Reproducing published S&L numbers requires S&L's simplified atmosphere (`adc` routine) and g = 32.17 ft/s²; matching JSBSim requires US 1976 atmosphere and JSBSim's gravity. Both are needed and selected per scenario.
7. **Quaternion linearization pitfall.** Linearizing in quaternion coordinates yields a redundant state and a spurious zero eigenvalue. Linearize in Euler-angle flight variables instead.
8. **Missing from the brief:** unit policy; error-metric definitions (attitude error must handle q ≡ −q); tolerances derived from error analysis; termination at ground contact; NaN/Inf guards; definition of reproducibility (bitwise vs tolerance); data provenance; sanitizer CI; controller update rate vs integrator step.

---

## 3. Goals, non-goals, anti-goals

### 3.1 Goals
- G1. Correct, tested 6-DOF rigid-body dynamics with explicit, tested frame conventions.
- G2. Integrators: explicit Euler, RK4, symplectic (semi-implicit) Euler, Lie–Euler and RKMK4 (exp-map attitude). Measured convergence orders and work–precision comparison.
- G3. Validation suite of analytic, conservation, and balance tests, run in CI.
- G4. F-16 nonlinear model reproducing the S&L published state-derivative check case.
- G5. Trim and linearization reproducing published S&L trim and linear-model results; modal analysis.
- G6. Matched-model JSBSim cross-validation with a quantified discrepancy budget.
- G7. Altitude + airspeed autopilot with actuator dynamics/limits and anti-windup; disturbance rejection under wind and discrete gusts; metrics reported.
- G8. One-command reproducible demo producing all plots and the numbers in `docs/results.md`.
- G9. Benchmarks of step cost and real-time factor.

### 3.2 Non-goals
- 3D graphics, GUI, real-time pacing, hardware-in-the-loop, networking.
- Landing gear, ground handling, takeoff/landing. Simulations terminate at h ≤ 0.
- Rotating/ellipsoidal Earth in our model; navigation over long distances.
- Post-stall, spin, or any flight outside the S&L table domain (α ∈ [−10°, 45°]; β limited by the narrowest table, verify on transcription; subsonic per data validity).
- Sensors, noise, state estimation (EKF), guidance.
- Flexible structures, multibody, fuel burn / time-varying mass.
- Dryden/von Kármán continuous turbulence (stretch only).
- Heading/lateral autopilot, LQR, coordinated-turn trim (stretch only).

### 3.3 Anti-goals (failure even if it "works")
- Plots that look plausible but have no quantitative check.
- Any tuning factor added to make results match JSBSim or S&L.
- Tolerances loosened without a written justification in the test and in `docs/design.md`.
- Claims of fidelity beyond what a listed validation supports.
- Non-reproducible results (unseeded randomness, results not regenerable from committed configs).

---

## 4. Conventions (normative)

Every item in this section has at least one dedicated test (§10.2).

### 4.1 Frames
| Frame | Symbol | Definition |
|---|---|---|
| World / inertial | N | Flat-Earth NED: origin fixed at ground level, x North, y East, z Down. Treated as inertial. Altitude h = −p_D. |
| Body | B | Origin at CG. x forward (nose), y right wing, z down (FRD). |
| Wind (air-relative) | W | x along air-relative velocity. α, β defined below. |
| JSBSim structural | S | x aft, y right, z up, inches. **Only** inside the JSBSim XML generator. |

- Direction cosine matrix `C_NB` maps body components to NED components: `v_N = C_NB v_B`; `C_BN = C_NBᵀ`.
- Euler angles: 3-2-1 sequence (ψ yaw, θ pitch, φ roll), `C_NB = R_z(ψ) R_y(θ) R_x(φ)`. Used **only** for I/O, initial conditions, linearization, and controller feedback; never integrated.
- Air-relative velocity: `v_air,B = C_BN (v_N − w_N)` with components (u, v, w).
  - `V = ‖v_air,B‖`, `α = atan2(w, u)`, `β = asin(v / V)` (guarded for V → 0).
- Body rates `ω_B = (p, q, r)`: angular velocity of B relative to N, in B components.
- Positive control deflections follow S&L: trailing-edge-down elevator positive (produces nose-down moment), etc. Exact sign table recorded in `docs/design.md` when transcribing Appendix A and tested by sign tests (§10.2).

### 4.2 Quaternions
- `q = q_NB = [q0, q1, q2, q3]` (scalar first), Hamilton product, unit norm.
- Vector rotation: `[0, v_N] = q ⊗ [0, v_B] ⊗ q*`.
- Kinematics: `q̇ = ½ q ⊗ [0, ω_B]`.
- Eigen pitfall: `Eigen::Quaterniond(w, x, y, z)` constructor is scalar-first but `.coeffs()` storage is `(x, y, z, w)`. State vectors store scalar-first explicitly; conversion functions are the only place Eigen quaternions are constructed, and are tested.
- Sign ambiguity: q and −q are the same rotation. All attitude comparisons use the error angle
  `θ_err = 2·atan2(‖vec(q_refᶜ ⊗ q)‖, |scalar(q_refᶜ ⊗ q)|)` (robust near zero, sign-invariant).
- Normalization policy: for integrators that integrate q̇ in ℝ⁴ (Euler, RK4), q is renormalized after each **full** step (projection). The norm drift *without* renormalization is measured and reported. Lie-group integrators preserve the norm by construction and do not renormalize.
- Euler-angle extraction clamps `asin` argument to [−1, 1]; behavior at θ = ±90° is defined and tested (ψ, φ split is convention: set φ = 0).

### 4.3 Units and constants
- SI internally (m, s, kg, N, rad). Imperial ↔ SI conversion constants defined once (`ft = 0.3048 m` exact, `slug = 14.59390294 kg`, `lbf = 4.4482216152605 N`).
- Gravity is a model parameter: constant `g0`, or inverse-square `g(h) = GM/(R+h)²`. Defaults: `g0 = 9.80665 m/s²`. S&L validation scenarios use `g0 = 32.17 ft/s²` (S&L value). JSBSim comparison uses inverse-square with JSBSim planet constants.
- Angles internally in radians; degrees only in configs/plots where labeled `_deg`.

### 4.4 State vectors
Rigid body, N = 13:
```
x_rb = [ p_N (3) | v_N (3) | q_NB (4) | ω_B (3) ]
```
F-16, N = 17:
```
x_ac = [ x_rb (13) | P (engine power level, %) | δe, δa, δr (actuator positions, rad) ]
```
Optional budget augmentation (for balance tests, §8.4): `+ [W_nc (1) | J_N (3)]` = work done by non-gravitational forces and angular impulse in N.

Inputs (F-16): `u = [throttle ∈ [0,1], δe_cmd, δa_cmd, δr_cmd]`. In S&L-reproduction mode actuators are bypassed (surface position = command) to match the book's model.

---

## 5. Mathematical model

### 5.1 Rigid-body equations of motion
```
ṗ_N = v_N
v̇_N = (1/m) C_NB F_B + g_N(h),           g_N = [0, 0, g(h)]
q̇   = ½ q ⊗ [0, ω_B]
ω̇_B = I⁻¹ ( M_B − ω_B × (I ω_B + h_rotor,B) )
```
- `F_B`, `M_B`: sum of non-gravitational body-frame forces/moments about the CG.
- `I`: full 3×3 inertia tensor about CG in B (products of inertia allowed); `I⁻¹` precomputed.
- `h_rotor,B`: constant rotor angular momentum (F-16 engine, S&L `hx`, along body x). Zero for generic bodies.
- Mass and inertia are constant.

### 5.2 Force/moment composition (generic engine)
A scenario's forces are a list of callables `ForceMoment(double t, const State&) → {F_B, M_B}` (or NED-expressed forces converted at the boundary). Provided generic models (used by validation scenarios):
- Uniform gravity (built in), constant body force/torque, linear drag `F = −c v`, quadratic drag `F = −k ‖v‖ v`, linear spring `F = −K p` (isotropic or diagonal anisotropic).

### 5.3 F-16 model (S&L Appendix A)
The S&L low-fidelity tabulated model (derived from NASA TP-1538). **Transcription sources** (owner has no copy of the book):
- **Primary transcription:** the numeric tables and constants in an open reproduction of the S&L model, e.g. AeroBenchVV (Heidlauf et al., MATLAB, `github.com/pheidlauf/aerobenchvv`) or its Python port (`github.com/stanleybak/AeroBenchVVPython`).
- **Independent cross-check:** a second reproduction that did not derive from the first (e.g. F16Model.jl `github.com/isrlab/F16Model.jl`, or the Archimedes "Subsonic F-16" tutorial). Every table and constant is diffed programmatically between the two; disagreements are resolved against NASA TP-1538 (NTRS 19800005879) where it covers the value, and documented.
- Only numbers are transcribed; no code is copied (several sources are GPL). Each table in the TOML cites TP-1538 and the reproduction(s) it was checked against.

Model content:
- **Geometry and mass** (verify every value on transcription and record source): W = 20,500 lbf; Ixx = 9,496, Iyy = 55,814, Izz = 63,100, Ixz = 982 slug·ft²; S = 300 ft²; b = 30 ft; c̄ = 11.32 ft; x_cg,ref = 0.35 c̄; h_x = 160 slug·ft²/s.
- **Aerodynamics** (body axes): tabulated `Cx(α, δe)`, `Cy(β, δa, δr)`, `Cz(α, β, δe)`, `Cl`, `Cm(α, δe)`, `Cn`, control-effectiveness tables, and damping-derivative tables vs α (CXq, CYr, CYp, CZq, Clr, Clp, Cmq, Cnr, Cnp) with nondimensional rates `p b/(2V)`, `q c̄/(2V)`, `r b/(2V)`. CG shift: `Cm += Cz (x_cg,ref − x_cg)`, `Cn −= Cy (x_cg,ref − x_cg) c̄/b` (as in S&L).
- **Interpolation:** bilinear/linear with the **same edge behavior as the S&L reference code** (linear extrapolation one cell beyond the table). A runtime flag records any lookup outside the table domain; validation and comparison runs assert none occurred. A **breakpoint monitor** (optional, test/diagnostic builds of the runner) records which table cells each lookup used, so tests can assert that no breakpoint was crossed (used by LN2 and convergence diagnostics).
- **Table domains** are whatever Appendix A tabulates (e.g. control-effectiveness tables may cover a narrower β range than the static tables); the model's valid envelope is the intersection, recorded in `docs/design.md`.
- **Engine:** throttle → commanded power via S&L "throttle gearing"; first-order power lag with S&L's power-dependent time constant; thrust from idle/military/maximum tables vs altitude and Mach; thrust along body x through CG; engine angular momentum `h_x` in the moment equation.
- **Atmosphere:** two implementations selectable per scenario:
  - `sl_adc`: S&L's simplified troposphere/stratosphere routine (used to reproduce book results).
  - `us76`: U.S. Standard Atmosphere 1976, 0–20 km (used for JSBSim comparison and general scenarios).
- **Actuators** (control milestone): first-order lag + rate limit + position limit per surface (S&L values: elevator ±25°, aileron ±21.5°, rudder ±30°; lag ≈ 20.2 rad/s; rate limits per S&L autopilot chapter — verify on transcription). Throttle limited to [0, 1].
- **Data source of truth:** `data/f16/f16_sl.toml` containing every table with units, breakpoints, and a citation (book section/page and TP-1538 reference). C++ loads it; the JSBSim XML generator reads the same file.

### 5.4 Wind and disturbances
- Constant wind vector `w_N`.
- Discrete gust: "1-cosine" shape (MIL-F-8785C), amplitude `V_m`, length `d_m` along a chosen NED or body axis, onset time **aligned to an integration step boundary**.
- Wind enters only through air-relative velocity. Spatial wind gradients (rotational gust effects) are neglected and documented.
- Dryden turbulence: stretch (seeded RNG if implemented).

### 5.5 Termination and guards
- Stop on: final time; `h ≤ 0` (ground contact); any non-finite state (NaN/Inf), reported with time and offending state; out-of-domain aero lookup (configurable: warn or stop).

---

## 6. Numerical integration

### 6.1 Integrators
All integrators are free function templates over `Eigen::Matrix<double, N, 1>` and a derivative callable `f(t, x, u)`. Inputs `u` are held constant over each step (zero-order hold).

| Name | Order | Notes |
|---|---|---|
| Explicit Euler | 1 | Baseline. |
| RK4 (classical) | 4 | Quaternion renormalized after each full step. |
| Symplectic (semi-implicit) Euler | 1 | `v_{k+1} = v_k + h a(p_k)`, `p_{k+1} = p_k + h v_{k+1}`. Used on separable translational problems (harmonic oscillator) to show bounded energy error. Not applied to rotation (not symplectic there); documented. |
| Lie–Euler | 1 | Attitude: `q_{k+1} = q_k ⊗ exp(½ h ω_k)`; other states Euler. |
| RKMK4 (Runge–Kutta–Munthe-Kaas) | 4 | Attitude parametrized as `q(t) = q_k ⊗ exp(½ Θ(t))`, Θ ∈ so(3) ≅ ℝ³ in body coordinates, Θ(0) = 0. Because the update multiplies on the right with body rates, the stage equation is `Θ̇ = dexp⁻¹_{−Θ}(ω_B) ≈ ω_B + ½ Θ×ω_B + (1/12) Θ×(Θ×ω_B)` (the right-Jacobian inverse; note the **+½** sign — the opposite sign silently degrades the method to order 2 on R2/R3 while R1 still passes). Classical RK4 tableau applied to Θ and to the vector-space states; step ends with `q_{k+1} = q_k ⊗ exp(½ Θ_{k+1})`. Norm preserved to roundoff without renormalization. |

The quaternion exponential uses a series expansion for small angles (‖θ‖ < 1e-4) to avoid 0/0.

### 6.2 Multi-rate structure
- Integrator step `dt` (default 0.01 s for aircraft; varies in studies).
- Controller runs at fixed rate (default 50 Hz); `dt` must divide the controller period exactly (validated at scenario load).
- Discontinuous events (gust onset, step inputs) occur only at step boundaries (validated at scenario load).

### 6.3 Convergence methodology (normative)
- Fixed final time T. Step sequence `h_k = h_0 / 2^k`, k = 0…K (K ≥ 7).
- Error at T against reference: analytic where available; otherwise a fine-step **RK4 (or RKMK4)** solution at `h_ref = h_K / 16`, regardless of the method under test. The reference's own error (estimated from its self-convergence between `h_K/8` and `h_K/16`) must be < 1e-3 × the finest-level error of the method under test, otherwise the reference is refined further. (A same-method reference biases first-order methods' observed order upward, e.g. ≈ 1.1 at the finest level.)
- Error metrics: position/velocity/ω Euclidean norm; attitude via `θ_err` (§4.2); energy and momentum relative drifts.
- Observed order: `p_k = log2(e(h_k)/e(h_{k+1}))`. Report the least-squares slope over the **asymptotic window**: levels whose error is > 1e3·ε_machine·scale (above roundoff floor) and whose successive `p_k` agree within 10%.
- Every convergence study outputs: error-vs-h log-log plot with reference slopes, observed-order table, and **work–precision** plot (error vs. number of derivative evaluations).
- Aircraft-model convergence is measured and reported; reduced order due to non-smooth tables is expected and explained, not asserted as 4.

---

## 7. Aircraft analysis

### 7.1 Trim
- Steady wings-level flight at specified (V, h, γ, x_cg): unknowns `z = [α, throttle, δe]`; β = δa = δr = φ = p = q = r = 0, θ = α + γ, engine power at its commanded steady state.
- Residual: `r(z) = [V̇, α̇, q̇]` evaluated through the nonlinear model (scaled to comparable magnitudes).
- Solver: Levenberg–Marquardt / damped Newton with central-difference Jacobian and backtracking. Converged when `‖r‖∞ < 1e-10` (scaled) within 50 iterations; failure reported with last iterate.
- Bounds respected: throttle ∈ [0,1], δe within limits, α within table domain; out-of-bounds solution is reported as "no trim".
- Should (not required): general trim with β, δa, δr for coordinated turns (S&L turn-coordination constraint).

### 7.2 Linearization
- Flight-variable chart (S&L state order): `x13 = [V, α, β, φ, θ, ψ, p, q, r, p_N, p_E, h, P]`, inputs `[throttle, δe, δa, δr]`.
- `f13(x13, u)`: map x13 → model state, evaluate model derivative, map derivative back (body acceleration → V̇, α̇, β̇; body rates → Euler rates via the kinematic matrix). This same mapping is used to compare against the S&L check case (§8.5), so it is validated independently.
- Jacobians `A = ∂f13/∂x`, `B = ∂f13/∂u` by central differences with per-variable step `δ_i = ε·max(|x_i|, s_i)`, ε ≈ 6e-6 (≈ ε_machine^(1/3)), `s_i` a per-variable scale.
- Outputs: A, B as JSON; longitudinal and lateral-directional sub-models; eigenvalues; mode identification (short period, phugoid, Dutch roll, roll subsidence, spiral) with ωn, ζ, time constants.

### 7.3 Control (required: altitude + airspeed hold)
- Discrete PID, sampled at controller rate, with: derivative on measurement, first-order derivative filter, output saturation, and **back-calculation anti-windup** (tracking time constant `T_t`). Anti-windup can be disabled to demonstrate windup.
- Longitudinal autopilot (cascade): `h_cmd → [PI] → θ_cmd (limited ±15°) → [P] → q_cmd → [PI] → δe_cmd`. Airspeed: `V_cmd → [PI] → throttle` with throttle limits and engine lag.
- Design workflow: gains designed in Python on the linear model (root locus/Bode), loop-at-a-time gain and phase margins reported, then verified on the nonlinear model.
- Stretch: heading hold (ψ → φ_cmd → δa, rudder coordination), LQR comparison.

---

## 8. Validation strategy

### 8.1 Validation matrix — what each check proves

| Layer | Check | Independent reference | Proves | Cannot prove |
|---|---|---|---|---|
| Kinematics | Rotating body-fixed force (T5), constant ω (R1), symmetric top (R2) | Closed form | Frame-transform direction, quaternion kinematics, integrator order | Inertia coupling |
| Rigid body | Asymmetric torque-free (R3), H_N & T conservation | Jacobi elliptic closed form; conservation laws | Euler's equations incl. gyroscopic term; coupling of ω and attitude | Force models |
| Integrators | Convergence studies, work–precision | Analytic solutions | Observed order, cost efficiency | Model correctness |
| Linearizer | Spin-stability eigenvalues (L1) | Closed form | Jacobian machinery | Aircraft specifics |
| Aircraft EOM | Work–energy (B1), angular-impulse balance (B2), wind Galilean invariance (W1), mirror symmetry (S1) | Physical laws | Euler's equations and force/moment assembly (B1), attitude kinematics consistency (B2), wind handling (W1), lateral sign consistency (S1) | Aero data values; transform direction (owned by T5) |
| Aircraft data | Two-source table diff; S&L state-derivative check case (A1) | Independent reproductions; published numbers (via secondary sources) | Table transcription, model assembly | Integration over time |
| Trim / linear | S&L published trim and linear model (TR1, LN1, via secondary sources), linear-vs-nonlinear scaling (LN2) | Published numbers; Taylor theorem | Trim solver, linearization, mode identification | Dynamic response fidelity |
| Implementation | Matched JSBSim (J1–J3) | Independent simulator, same data | EOM, integration, atmosphere, table handling, conventions | Aero fidelity to real F-16 |
| Model form | Stock JSBSim F-16 appendix (J4) | Different data set | Scale of model-form uncertainty | Anything about correctness |

### 8.2 Translational cases (analytic)
| ID | Case | Reference | Role |
|---|---|---|---|
| T1 | Constant acceleration | Polynomial | Smoke test only (RK4 exact) |
| T2 | 3D linear drag + gravity | `v(t) = v∞ + (v0 − v∞) e^{−ct/m}`, p(t) integrated | Convergence (Euler, RK4) |
| T3 | Vertical quadratic drag, release from rest | `z(t) = (v_t²/g) ln cosh(g t / v_t)`, `v = v_t tanh(g t / v_t)` | Convergence; nonlinear force |
| T4 | Harmonic oscillator (anisotropic spring) | sin/cos; energy | Energy behavior: Euler grows, RK4 slow secular drift, symplectic bounded |
| T5 | Body spinning at constant ω about body z with constant body-x force | Closed-form rotating force integral | Catches `C_NB` vs `C_BN` error (wrong direction of rotation) |
| T0 | Zero force | v constant | Smoke test only |

### 8.3 Rotational cases
| ID | Case | Reference |
|---|---|---|
| R1 | Torque-free spin about a principal axis | `q(t) = q0 ⊗ exp(½ ω t)` (RK4 not exact → real convergence test) |
| R2 | Torque-free axisymmetric body (I1 = I2 ≠ I3) | `q(t) = exp(½ Ω_p t ĥ_N) ⊗ q0 ⊗ exp(−½ λ t e3)`, `Ω_p = ‖L‖/I1`, `λ = ω3 (I3 − I1)/I1`. Implementation must verify the reference itself: numerically differentiate the formula and confirm it satisfies `q̇ = ½ q ⊗ ω` to 1e-9. |
| R3 | Torque-free asymmetric body (I1 < I2 < I3, L² > 2E I2), ω2(0) = 0 | ω(t) via Jacobi elliptic functions (Landau & Lifshitz §37). Reference self-check: satisfies Euler's equations (finite-difference residual < 1e-9) and conserves E, ‖L‖. Attitude checked via conservation of `H_N = C_NB I ω_B` (vector) and fine-step self-convergence. |
| R4 | Intermediate-axis (Dzhanibekov) flip | Demonstration plot + conservation; with L1 below gives the stability explanation |
| R5 | Torque-free body with rotor `h_rotor` | Conservation of total `H_N = C_NB (I ω_B + h_rotor)` |

Jacobi elliptic functions sn, cn, dn: implemented with the AGM/descending Landen method (Abramowitz & Stegun 16.4); validated against committed reference values generated by `scipy.special.ellipj` and the identities `sn² + cn² = 1`, `dn² + k² sn² = 1`.

Conservation metrics reported for every rotational case and integrator: max relative energy drift, max relative ‖H_N(t) − H_N(0)‖, max |‖q‖ − 1| (with and without renormalization for RK4/Euler).

### 8.4 Aircraft consistency laws
| ID | Check | Criterion |
|---|---|---|
| B1 | Work–energy including rotation: `E = ½m‖v_N‖² + ½ω_BᵀIω_B + m·Φ(h)` vs `W_nc = ∫ [(C_NB F_B)·v_N + M_B·ω_B] dt` (augmented state). Gyroscopic terms do no work (`ω·(ω×(Iω+h)) = 0`), so this exercises Euler's equations and the force/moment assembly. It does **not** test the direction of `C_NB` (both sides use it) — T5 and B2 own that. | `|ΔE − W_nc| / max|ΔE|` converges at the integrator's observed order and is < 1e-6 at dt = 0.01 s for a 60 s maneuvering run |
| B2 | Angular impulse: `ΔH_N` vs `∫ C_NB M_B dt` (augmented state; includes rotor) | Same criterion |
| W1 | Galilean invariance: trimmed flight in a constant uniform wind | Air-relative states identical (≤ 1e-9 rel.) to the no-wind run; ground track offset by `w_N t` |
| S1 | Mirror symmetry with `h_x = 0`: ±δa doublets, run on a **test-only symmetrized copy** of the aileron/rudder effectiveness tables (made even in β). The published tables are not exactly symmetric in β, so the production model is not expected to mirror exactly. | Symmetrized data: lateral states exact negatives, longitudinal identical (≤ 1e-12 rel.). Production data and `h_x ≠ 0`: asymmetry measured and reported, not asserted |
| L1 | Linearizer on torque-free spin about principal axis 1 at rate Ω | Eigenvalues match `λ² = −Ω² (I1 − I2)(I1 − I3)/(I2 I3)` to 1e-6 rel. (stable for major/minor axis, unstable for intermediate) |

### 8.5 Published-data checks (F-16)
| ID | Check | Criterion |
|---|---|---|
| A1 | S&L state-derivative check case (book Eq./Table 3.5-2: V = 500 ft/s, α = 0.5 rad, β = −0.2 rad, x_cg = 0.4, P = 90 %, plus the book's remaining state/control values) | Every derivative component within 1e-3 relative (book prints ~4–6 significant figures). Any component differing by > 1e-4 relative is investigated and explained in `docs/results.md`. |
| TR1 | S&L published trim results (level flight, sea level and altitude cases as tabulated in the book; x_cg = 0.35) | α, throttle, δe match to the book's printed precision (≤ 1 in the last printed digit, or ≤ 0.5 % rel., whichever is larger) |
| TR2 | Stationarity: simulate 10 s from trim with dt = 0.01 s | |ΔV| < 1e-3 ft/s, |Δα|, |Δθ| < 1e-6 rad, |Δh| < 1e-2 ft (tolerances to be confirmed against the unstable/neutral mode growth from 1e-10 residuals; any revision documented) |
| LN1 | S&L published linear model (502 ft/s, sea level, x_cg = 0.35 example) | A, B entries within 1 % or 1 in the last printed digit; eigenvalues within 1 % |
| LN2 | Linear vs nonlinear elevator step of δ and δ/10, 5 s. The trim δe sits near the δe = 0 table breakpoint, so the step is **trailing-edge-up (more negative δe)**, δ = 0.5°, and the run asserts via the interpolation breakpoint monitor (§5.3) that no α/β/δ breakpoint cell is crossed (crossing a kink makes the error scale ~linearly, not quadratically) | Max-error ratio e(δ)/e(δ/10) in [70, 130] |
| LN3 | Finite-difference step sweep (ε from 1e-2 to 1e-10) | Plot shows truncation/roundoff "V"; chosen ε lies in the flat region |

**Reference-value sourcing (no book access).** A1, TR1 and LN1 reference numbers are the S&L published values **as reproduced in at least two independent secondary sources** (open reproductions' test files/tutorials, papers that tabulate them). Each value in `docs/results.md` lists its sources. If only one source is found, the check still runs but is labeled "single-source"; if sources disagree, the disagreement is reported and the check is evaluated against each. `docs/design.md` states plainly that these are secondary reproductions of book values, which is weaker independence than direct transcription. If the owner later gains access to the book (e.g. via a library), these become direct transcriptions with page citations; no other change is needed. Values in this PRD that were not verified against a source are marked "verify".

### 8.6 JSBSim cross-validation (matched model)
**Principle:** remove every known difference except the implementations, then quantify what remains.

Generated JSBSim aircraft (`python/jsbsim_compare/gen_xml.py` → `reference/jsbsim/aircraft/f16_sl/f16_sl.xml`), from `data/f16/f16_sl.toml`:
- Aero tables as JSBSim `<function>`/`<table>` in body axes X/Y/Z and moment axes, using `aero/qbar-psf`, `aero/bi2vel`, `aero/ci2vel`, `aero/alpha-rad`, `aero/beta-rad`.
- **Table edges:** JSBSim clamps at table edges whereas S&L extrapolates one cell → comparison runs are constrained to stay inside the table domain (asserted on both sides).
- **CG/AERORP:** AERORP at `x_cg,ref`, CG at `x_cg` in the structural frame (x aft, inches) so JSBSim's moment transfer reproduces S&L's Cm/Cn CG shift. Verified by A1-style derivative comparison at t = 0.
- **Inertia:** `<mass_balance negated_crossproduct_inertia="false">` (JSBSim's default negates products of inertia). Ixz sign verified by comparing initial ṙ response to a roll-rate initial condition.
- **Engine:** thrust implemented as a body-x external force computed from the same thrust tables at constant power P (trim power). Engine power-lag dynamics are **excluded** from comparison cases (constant throttle); documented.
- **Rotor gyroscopic term:** `h_x` set to 0 in both simulators for comparison runs (JSBSim has no direct equivalent); documented.
- **Planet:** custom planet file with `rotation_rate = 0`, spherical, `simulation/gravity-model = 0` (standard point-mass). Our model runs with matching inverse-square gravity. Remaining Earth-curvature effect estimated analytically and reported.
- **Atmosphere:** US 1976 on both sides.
- **Integration:** JSBSim run at two rates (e.g. 1/1200 s and 1/2400 s) to show its discretization error is ≪ the observed discrepancy; our side runs RK4 at a step where its own error is likewise negligible. JSBSim integrator settings recorded in metadata.

Cases:
| ID | Case | Criterion |
|---|---|---|
| J1 | Point check: JSBSim initialized at our trim state; report JSBSim's initial accelerations | `|u̇|,|ẇ|` < 1e-3 ft/s², `|q̇|` < 1e-5 rad/s² (initial targets; exceeding them blocks progress until root-caused) |
| J2 | Each simulator's own trim at the same (V, h) | α, δe, throttle agree within 0.1 % (or explained) |
| J3 | Open-loop responses from each simulator's own trim, 20 s: elevator doublet, aileron pulse, rudder doublet | Compare deviations from trim (ΔV, Δα, Δq, Δθ, Δh, Δβ, Δp, Δr, Δφ). Target: max difference ≤ 2 % of each signal's peak deviation |
| J4 (appendix) | JSBSim stock `f16` at a matched condition | Trim values only: the M4.0 spike showed the stock FCS overwrites direct surface commands, so open-loop responses would require editing a copy of the GPL model (out of scope). Reported as **model-form** difference; no pass/fail beyond "runs and is explained" |

**Amendment from the M4.0 spike (2026-11-08, `spikes/jsbsim/FINDINGS.md`, jsbsim 1.3.3):** the primary approach is feasible (custom planet via `fdm.load_planet`, `simulation/gravity-model = 0`, `negated_crossproduct_inertia="false"`, direct surface drive, JSBSim trim). JSBSim has no RK4; its integrators are set explicitly (recorded in metadata) and its own discretization error is bounded by comparing two small steps (e.g. 1/2400 s and 1/4800 s). Matched-model trim requires a dummy zero-power engine so FGTrim can vary throttle, with thrust supplied by the external force.

**Discrepancy budget** (required section in `docs/results.md`): each candidate source (gravity model, curvature, atmosphere implementation, table edge handling, integrator, sampling alignment, residual trim error) is toggled or bounded individually and its contribution to J3 differences quantified. Remaining unexplained difference must be stated plainly.

### 8.7 Control performance
Scenarios at the S&L 502 ft/s sea-level-ish trim, with actuators and engine lag enabled. Initial targets (finalized at the M7 review; changes documented):

| ID | Scenario | Metric targets |
|---|---|---|
| C1 | +100 ft altitude step | Overshoot ≤ 10 %; 2 % settling ≤ 25 s; steady-state |e| ≤ 1 ft; airspeed deviation ≤ 5 ft/s |
| C2 | +1,000 ft altitude step (saturating θ_cmd and throttle) | With anti-windup: overshoot ≤ 15 %; without: windup visibly worse (plotted side-by-side) |
| C3 | +50 ft/s airspeed step at constant altitude | Overshoot ≤ 10 %, altitude deviation ≤ 20 ft |
| C4 | Vertical 1-cos gust (amplitude 10 ft/s) and constant 20-kt headwind change | Peak altitude deviation and recovery time reported; returns within 2 % band |
| C5 | Linear-model loop margins | Gain margin ≥ 6 dB, phase margin ≥ 45° per loop (broken one loop at a time) |

Reported metrics: overshoot, 2 % settling time, steady-state error, peak actuator rate/position usage, linear-predicted vs nonlinear-measured step response overlay.

---

## 9. Architecture

### 9.1 Repository layout
```
CMakeLists.txt, CMakePresets.json
cmake/Dependencies.cmake        # FetchContent, pinned versions + hashes
src/flightsim/
  math/        quat.hpp, frames.hpp, interp.hpp, elliptic.hpp
  core/        state.hpp, rigid_body.hpp/.cpp, integrators.hpp, forces.hpp
  env/         atmosphere.hpp/.cpp (sl_adc, us76), wind.hpp/.cpp
  aircraft/    f16_data.hpp/.cpp (TOML loader), f16_model.hpp/.cpp (aero, engine, actuators)
  analysis/    trim.hpp/.cpp, linearize.hpp/.cpp, flight_vars.hpp/.cpp (x13 chart)
  control/     pid.hpp/.cpp, autopilot.hpp/.cpp
  sim/         scenario.hpp/.cpp (TOML), runner.hpp/.cpp, output.hpp/.cpp (CSV + JSON metadata)
apps/flightsim.cpp              # CLI
tests/                          # GoogleTest; one file per module + validation/
bench/                          # Google Benchmark
data/f16/f16_sl.toml            # single source of truth for F-16 data, with citations
scenarios/*.toml                # every scenario used in docs/results.md
python/                         # pyproject.toml, pinned requirements
  flightsim_analysis/           # plots, convergence fits, report generation
  jsbsim_compare/               # gen_xml.py, run_jsbsim.py, compare.py
  control_design/               # linear-model design, margins
reference/jsbsim/               # generated aircraft/planet XML + committed reference CSVs
scripts/demo.sh                 # one-command demo
docs/PRD.md, docs/design.md, docs/results.md
```
Modules are folders in one static library (`flightsim_core`), not separate libraries. No abstract base classes are required; force models and derivative functions are callables; integrators are templates.

### 9.2 Key interfaces (shape, not final code)
```cpp
template <int N> using Vec = Eigen::Matrix<double, N, 1>;

// Integrators: free functions, ZOH input u.
template <int N, class F, class U> Vec<N> euler_step(F&& f, double t, const Vec<N>& x, const U& u, double h);
template <int N, class F, class U> Vec<N> rk4_step  (F&& f, double t, const Vec<N>& x, const U& u, double h);
// Lie-group variants know the attitude slice layout of the state.

struct ForceMoment { Eigen::Vector3d F_B, M_B; };
using ForceModel = std::function<ForceMoment(double t, const RigidBodyState&)>;

struct F16Model {             // pure: no hidden state
  F16Data data; AtmosphereKind atm; GravityModel grav; WindModel wind; bool actuators;
  Vec<17> deriv(double t, const Vec<17>& x, const F16Input& u) const;
};
```
- Derivative functions are pure (no mutable state), so RK stages, Jacobians, and trim can call them freely.
- The runner owns time, controller scheduling, event alignment, termination, and output.
- No heap allocation inside the step loop (verified by a test that counts allocations via a replaced `operator new` in a dedicated test binary).

### 9.3 CLI
```
flightsim run <scenario.toml> --out <dir>          # trajectory CSV + metadata JSON
flightsim trim <scenario.toml> --out <file.json>
flightsim linearize <scenario.toml> --out <file.json>
flightsim converge <study.toml> --out <dir>         # convergence study tables
```
- CSV: one row per output sample; column names carry units (`t_s, pN_m, ..., alpha_rad, V_mps`), includes derived quantities (Euler angles, V, α, β, h, Mach, q̄, E, H_N, ‖q‖).
- Metadata JSON: git SHA, build type, compiler, scenario file hash, integrator, dt, all resolved parameters.

### 9.4 Scenario file (illustrative)
```toml
# Torque-free asymmetric body: validates Euler's equations against Jacobi elliptic solution.
[sim]
integrator = "rk4"        # euler | rk4 | symplectic_euler | lie_euler | rkmk4
dt_s = 0.001
t_final_s = 20.0
output_every = 10

[body]
mass_kg = 1.0
inertia_kgm2 = [[1.0,0,0],[0,2.0,0],[0,0,3.0]]

[initial]
position_m = [0,0,0]
velocity_mps = [0,0,0]
euler_deg = [0,0,0]
omega_radps = [0.5, 0.0, 1.0]   # L² = 9.25 > 2E·I2 = 6.5 → Landau–Lifshitz branch used by R3; ω2(0) = 0

[[forces]]            # none for torque-free
```

### 9.5 Dependencies (all pinned; exact versions chosen in M0 and recorded)
- C++: Eigen 3.4.x, GoogleTest 1.15.x, Google Benchmark 1.9.x, toml++ 3.4.x — via `FetchContent` with `URL_HASH`.
- Python ≥ 3.11: numpy, scipy, plotly, kaleido (static PNG export), python-control, jsbsim (pinned), pytest. Locked in `python/requirements.lock`.
- Compilers: GCC ≥ 12, Clang ≥ 15, AppleClang ≥ 15. Warnings `-Wall -Wextra -Wpedantic -Wconversion -Werror` (or equivalent).

### 9.6 Reproducibility
- Same binary + same scenario ⇒ bitwise-identical CSV (tested). No threads in the simulation path; any RNG is seeded from the scenario.
- Across compilers/platforms, results agree within test tolerances, not bitwise (FMA contraction and libm differences). `-ffp-contract=off` set for determinism within a platform.
- Every number and plot in `docs/results.md` is regenerated by `scripts/demo.sh` from committed scenarios; `docs/results.md` tables are produced by the report generator from a metrics JSON, not typed by hand.

---

## 10. Testing requirements

### 10.1 Structure
- GoogleTest, registered with CTest, labels: `unit` (fast, < 10 s total in Debug), `validation` (analytic/convergence/published-data, < 2 min Release), `jsbsim` (requires Python + jsbsim; separate CI job).
- Python: pytest for the XML generator (round-trip: generated XML table values equal TOML values), convergence fitting, and report generation smoke test.

### 10.2 Required unit tests (non-exhaustive)
- Quaternion: Hamilton product vs hand-computed values; `q ⊗ q* = 1`; rotation of basis vectors by 90° about each axis; `C_NB(q)` orthonormal with det = +1; Euler↔quaternion↔DCM round trips over a grid (excluding |θ| > 89.9°) and defined behavior at gimbal lock; error angle invariant under q → −q; exp/log round trip incl. small-angle branch.
- Frames: a body-x unit vector at ψ = 90° points East; at θ = 90° points Up (−D); α > 0 for positive w; β > 0 for positive v.
- Sign tests (F-16): positive δe → negative pitching-moment increment; positive δa and δr signs per S&L convention; positive q → negative Cmq contribution (damping).
- Units: conversion constants; F-16 data loaded in SI equals TOML imperial × factor.
- Interpolation: exact at breakpoints; linear between; extrapolation behavior matches S&L code; out-of-domain flag set.
- Atmosphere: `us76` vs tabulated standard values at 0, 1, 5, 11, 15, 20 km (ρ, T, p, a) to 1e-4 rel.; `sl_adc` vs S&L routine values.
- Elliptic functions vs committed scipy values (1e-12 abs).
- PID: proportional/integral/derivative responses on scripted inputs; saturation; back-calculation tracks the saturated output; derivative-on-measurement has no kick on setpoint step.
- Scenario loader: rejects dt not dividing controller period, unaligned events, missing fields, unknown keys.
- Runner: NaN guard stops with diagnostic; ground-contact termination; determinism (two runs bitwise equal).

### 10.3 Validation tests (CI-asserted)
All cases in §8.2–§8.5 and J-series criteria are encoded as tests with the tolerances stated there. Convergence tests assert observed order (first-order methods ∈ [0.95, 1.05]; fourth-order ∈ [3.9, 4.1]) on the cases where the method is applicable and not exact:

| Method | Asserted on | Not asserted on (why) |
|---|---|---|
| Euler, RK4 | T2, T3, T4, R1, R2, R3 | — |
| Symplectic Euler | T4 | T2/T3 (velocity-dependent force, not separable), R* (not symplectic for rotation) |
| Lie–Euler, RKMK4 | R2, R3 | R1 (exact for constant ω: error at roundoff, no asymptotic window); T* (reduce to Euler/RK4) |

Additionally asserted: RK4/RKMK4 on R1 reproduce the exact solution to ≤ 1e-12 rad (RKMK4) and show order 4 (RK4).

### 10.4 Tolerance policy
Every tolerance carries a comment giving its origin: analytic error bound, observed error × safety factor (≤ 10×), or published precision. Changing a tolerance requires updating that comment and `docs/design.md`.

### 10.5 CI (GitHub Actions)
| Job | Content |
|---|---|
| build-test | ubuntu-24.04 {GCC, Clang} × {Debug, Release}; macOS AppleClang Release; `unit` + `validation` |
| sanitizers | Clang ASan + UBSan, `unit` + `validation` |
| format | clang-format check |
| python | pytest; generate plots from a short scenario set (smoke) |
| jsbsim | pinned jsbsim; regenerate J1–J3 references; compare to committed CSVs within tolerance |
| bench | build and run benchmarks once (results archived as artifact; not gated) |
| coverage (should) | llvm-cov report for `src/`; target ≥ 85 % line coverage on `core/`, `math/`, `analysis/` |

---

## 11. Performance benchmarks

Measured with Google Benchmark, Release, single thread; hardware recorded in `docs/results.md`. Targets are expectations, reported, not CI-gated.

| Benchmark | Target |
|---|---|
| Rigid-body RK4 step (gravity + quadratic drag) | ≤ 200 ns |
| F-16 derivative evaluation | ≤ 400 ns |
| F-16 RK4 step | ≤ 2 µs |
| 60 s closed-loop F-16 run, dt = 0.01 s, output disabled | ≤ 50 ms (≥ 1,000× real time) |
| Trim solve | ≤ 10 ms |
| Linearization (13×13 + 13×4) | ≤ 5 ms |
| Integrator comparison: cost per step for Euler / RK4 / RKMK4 | Reported, with work–precision plot |

---

## 12. Deliverables

1. C++20 library `flightsim_core` + CLI `flightsim`.
2. Integrators: Euler, RK4, symplectic Euler, Lie–Euler, RKMK4.
3. Test suite (unit + validation) green in CI on Linux and macOS, plus sanitizers.
4. Benchmarks with results table.
5. Committed scenarios for every result.
6. Python analysis: convergence plots, work–precision, conservation/balance plots, trim/linear/mode tables, JSBSim comparison plots, autopilot response plots (Plotly HTML + PNG).
7. F-16 data file with citations; generated JSBSim aircraft + planet files; committed JSBSim reference CSVs.
8. `scripts/demo.sh`: build → test → run all scenarios → generate plots → write `out/index.html` and regenerate `docs/results.md` tables. One command, from a clean clone, given the documented prerequisites.
9. `README.md`: one-paragraph pitch, 3–4 headline plots (convergence, Dzhanibekov + H_N conservation, JSBSim overlay, autopilot step), headline numbers, quickstart.
10. `docs/design.md`: architecture, conventions, derivations (EOM, quaternion kinematics, RKMK4, trim, linearization, autopilot structure), **decision log** (decision, alternatives, rationale), references with page numbers, known limitations.
11. `docs/results.md`: validation matrix with pass/fail and measured values, convergence tables, conservation drifts, published-data comparisons, JSBSim comparison + discrepancy budget, control metrics, benchmarks, J4 appendix.

---

## 13. Milestones

Each milestone ends with a review gate: all acceptance criteria met, CI green, docs for the milestone updated. Effort is in focused working days (estimates only).

### M0 — Foundations (1–2 d)
- CMake + presets, pinned FetchContent deps, CI skeleton (build-test, sanitizers, format), Python env lock, empty test passing, MIT license, README stub.
- **Accept:** clean clone builds and passes tests on all CI targets; `cmake --preset release && ctest --preset release` works locally on macOS.

### M1 — Rigid-body core + translational validation (3–5 d)
- Quaternion/frame utilities, 13-state EOM, force callables, Euler + RK4, runner, CSV/metadata output, TOML scenarios, CLI `run`.
- Tests: §10.2 quaternion/frame tests; T0–T5; determinism; NaN guard.
- **Accept:** T2/T3 observed orders within bounds; T5 passes; all conventions in §4 have a test.

### M2 — Rotational validation, convergence harness, extra integrators (5–7 d) → **MVP**
- R1–R5, elliptic functions, conservation metrics, convergence study runner (`flightsim converge`), symplectic Euler (T4), Lie–Euler, RKMK4, linearizer machinery + L1.
- **Accept:** all §10.3 order assertions pass for all five integrators; RK4 vs RKMK4 norm/conservation comparison table produced; L1 passes.

### M3 — Analysis tooling, demo, docs, benchmarks (3–4 d) → **Resume-ready core**
- Python plotting + report generator, `scripts/demo.sh`, benchmarks, `docs/design.md` (core sections), `docs/results.md` (core sections), README with headline plots.
- **Accept:** demo runs from clean clone in < 5 min (excluding dependency download) and regenerates every core result; zero-allocation test passes.

### M4 — F-16 model (5–8 d, **manual review: data transcription**)
- M4.0 **JSBSim feasibility spike (timeboxed 1 d, throwaway):** confirm via Python package: custom aircraft dir loads; custom planet with zero rotation loads (`LoadPlanet` or equivalent) ; `simulation/gravity-model` settable; integrator rates settable; surfaces drivable directly; stock `f16` surfaces drivable or not. Findings recorded in `docs/design.md`; §8.6 adjusted if needed.
- Data file with citations (primary transcription from one open reproduction, programmatically diffed against an independent second reproduction; disagreements resolved against TP-1538 and documented, §5.3); collect A1/TR1/LN1 reference values from ≥ 2 secondary sources; loader, interpolation, `sl_adc` + `us76` atmosphere, gravity models, engine, aero assembly, x13 flight-variable chart, termination on out-of-domain.
- Tests: A1, atmosphere, interpolation, sign tests, S1, W1, B1, B2.
- **Accept:** A1 within tolerance; B1/B2/S1/W1 pass; aircraft-model convergence measured and the reduced order explained in `docs/results.md`.

### M5 — Trim and linearization (5–7 d, **manual review: comparison with published values**)
- Trim solver, `flightsim trim`/`linearize`, mode identification, LN3 sweep.
- **Accept:** TR1, TR2, LN1, LN2 pass; mode table (ωn, ζ, τ) in results.

### M6 — JSBSim matched cross-validation (6–10 d, **highest risk; manual review**)
- XML + planet generator, JSBSim runner, comparison + plotting, discrepancy budget.
- **Accept:** J1–J3 criteria met or every exceedance root-caused and documented with quantified contribution; reference CSVs committed; jsbsim CI job green.

### M7 — Actuators and autopilot (6–10 d, **manual review: tuning**)
- Actuator models, PID with back-calculation anti-windup, cascaded altitude + airspeed hold, constant wind + 1-cos gust, linear-model design scripts (margins), nonlinear verification.
- **Accept:** C1–C5 met (or targets revised at review with written rationale); windup vs anti-windup comparison plot; linear-predicted vs nonlinear overlay.

### M8 — Stock JSBSim appendix and final polish (3–5 d)
- J4 appendix, full results regeneration, README final, decision log and limitations complete.
- **Accept:** a reader can trace every README claim to a test or a scenario + plot; fresh-clone demo succeeds on macOS and Linux.

### Stretch (only after M8)
Heading hold + coordinated-turn trim; LQR vs PID comparison; Dryden turbulence; smooth (C¹) table interpolation to show RK4 order recovery on the aircraft model.

---

## 14. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Table transcription errors | Wrong model silently | Two-source programmatic diff; A1 check case; TR1/LN1 as second net |
| No book access: reference values only via secondary sources, which may share errors | Weaker independence of A1/TR1/LN1 | ≥ 2 independent sources per value; label single-source checks; JSBSim matched comparison and physics-law checks (B1, B2, W1, S1) don't depend on book values; upgrade path if book obtained |
| JSBSim Python API lacks planet loading or direct surface control | J-series weaker | M4.0 spike before M6; fallback: quantify rotating-Earth/J2 effect analytically and by running our model with Coriolis-only perturbation estimate; compare deltas from own trims over short horizons |
| JSBSim inertia sign / structural-frame errors | Spurious discrepancies | Explicit `negated_crossproduct_inertia="false"`; J1 point check; Ixz-coupling sign test |
| Non-smooth tables degrade RK4 order on aircraft | "Doesn't converge at 4th order" | Predicted, measured, explained; stretch: C¹ interpolation |
| Trim fails to converge / multiple solutions | Blocks M5 | LM with bounds, good initial guess from published values, report failures explicitly |
| Unstable or lightly damped modes make TR2 drift | Flaky test | Short horizon; tolerances derived from measured residual growth |
| Autopilot tuning takes longer than estimated | Schedule | Design on linear model first; targets revisable at review with rationale |
| kaleido static export needs a Chrome install | Demo breaks on fresh machines | Document prerequisite; HTML outputs remain primary; PNGs committed for README |
| Book values printed to few digits | Ambiguous pass/fail | Tolerances tied to printed precision (§8.5) |
| Coding agent introduces sign/frame errors | Subtle wrong physics | Every convention has a test (§10.2); T5, S1, W1, B1/B2 catch classes of errors independent of data |

---

## 15. Open items to verify during implementation
1. Secondary sources for the S&L check case (A1), trim results (TR1), and 502 ft/s linear model (LN1): find ≥ 2 each; record which values are single-source.
2. S&L mass/geometry values listed in §5.3 (marked "verify") against both reproductions.
3. Actuator lag and rate limits (S&L autopilot chapter values as used in open reproductions, e.g. AeroBenchVV).
4. Licenses of the reproductions used (numbers only are transcribed; confirm nothing else is copied).
5. JSBSim Python exposure of `LoadPlanet` (added in JSBSim 1.2 via `FGFDMExec::LoadPlanet`), integrator property names, and whether stock `f16` surfaces can bypass its FCS (M4.0 spike).
6. Data provenance of JSBSim's stock `f16.xml` (read its header/references) for the J4 appendix narrative.

## 16. References
- Stevens, B. L., Lewis, F. L., Johnson, E. N. *Aircraft Control and Simulation*, 3rd ed., Wiley, 2016 — Appendix A F-16 model, Ch. 3 trim/linearization examples.
- Nguyen, L. T., et al. *Simulator Study of Stall/Post-Stall Characteristics of a Fighter Airplane with Relaxed Longitudinal Static Stability*, NASA TP-1538, 1979. https://ntrs.nasa.gov/citations/19800005879
- Heidlauf, P., Collins, A., Bolender, M., Bak, S. "Verification Challenges in F-16 Ground Collision Avoidance and Other Automated Maneuvers", ARCH 2018 — AeroBenchVV (S&L model reproduction), https://github.com/pheidlauf/aerobenchvv, https://github.com/stanleybak/AeroBenchVVPython
- F16Model.jl, https://github.com/isrlab/F16Model.jl; Archimedes "Subsonic F-16" tutorials, https://pinetreelabs.io/archimedes/tutorials/f16/
- Landau, L. D., Lifshitz, E. M. *Mechanics*, §37 (torque-free asymmetric top).
- Abramowitz, M., Stegun, I. *Handbook of Mathematical Functions*, §16.4 (Jacobi elliptic functions, AGM).
- Munthe-Kaas, H. "High order Runge–Kutta methods on manifolds", *Appl. Numer. Math.* 29 (1999); Iserles et al., "Lie-group methods", *Acta Numerica* (2000).
- Hairer, Lubich, Wanner. *Geometric Numerical Integration*, 2nd ed., Springer (symplectic Euler, energy behavior).
- U.S. Standard Atmosphere, 1976 (NOAA/NASA/USAF).
- MIL-F-8785C (1-cosine discrete gust).
- Åström, K. J., Murray, R. M. *Feedback Systems* (anti-windup by back-calculation).
- JSBSim reference manual and source (`FGInertial`, `FGMassBalance` `negated_crossproduct_inertia`), https://github.com/JSBSim-Team/jsbsim.
