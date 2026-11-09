# 6-DOF Flight Dynamics Simulator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Use superpowers:test-driven-development for every task: failing test first, then code.

**Goal:** Implement the approved PRD: a validated C++20 6-DOF rigid-body + F-16 flight dynamics simulator with convergence studies, trim/linearization, JSBSim cross-validation, and an altitude/airspeed autopilot.

**Architecture:** One static library `flightsim_core` (folders = modules), pure derivative functions over fixed-size Eigen vectors, integrators as free function templates, a thin CLI `flightsim`, GoogleTest suites split into `unit` and `validation` binaries, Python (uv-managed) for plots, reports, control design, and JSBSim.

**Tech Stack:** C++20, CMake ≥ 3.28 + presets, Eigen 3.4.0, GoogleTest 1.15.2, Google Benchmark 1.9.1, toml++ 3.4.0 (all FetchContent with URL_HASH); Python 3.12 via uv (numpy, scipy, plotly, kaleido, control, jsbsim, pytest).

**Spec:** `docs/PRD.md` (approved). Section numbers below (§x.y) refer to it. Executors read both.

## Global Constraints

- Namespace `flightsim`. Headers `.hpp`, sources `.cpp` under `src/flightsim/<module>/`. Include as `#include "flightsim/<module>/<file>.hpp"` (include root `src/`).
- SI units internally; radians internally. Imperial only at load/report boundaries; conversion constants only in `src/flightsim/math/units.hpp` (§4.3 exact values).
- Quaternion = `Eigen::Vector4d` scalar-first `[q0,q1,q2,q3]`, Hamilton, `q_NB` maps body→NED (§4.2). Never store `Eigen::Quaterniond` in state.
- Frames NED / FRD, `C_NB` maps body→NED, Euler 3-2-1 (§4.1).
- Rigid-body state layout `Vec<13> = [p_N(0..2) | v_N(3..5) | q_NB(6..9) | ω_B(10..12)]` (§4.4); index constants in `core/state.hpp`.
- Derivative functions are pure (no mutable state); no heap allocation in step loops.
- Warnings: `-Wall -Wextra -Wpedantic -Wconversion -Werror` on our targets; third-party via `SYSTEM`. `-ffp-contract=off` on all targets.
- Every tolerance carries a comment stating its origin (§10.4).
- Tests: GoogleTest. Binaries `unit_tests` (label `unit`) and `validation_tests` (label `validation`).
- Git: author `Eliot Boda <eliotboda@gmail.com>`. **No `Co-Authored-By` or "Generated with" lines** (owner's global rule). Logical commits, conventional prefixes (`feat:`, `test:`, `fix:`, `docs:`, `build:`, `ci:`).
- Do not copy third-party code (several F-16 sources are GPL). Numbers only, cited.

## Review Focus

1. Quaternion sign/order mistakes (`q⊗ω` vs `ω⊗q`, `C_NB` vs `C_BN`) that still pass symmetric tests — pinned by T5 (Task 1.3) and R2 (Task 2.2), which break on either error.
2. Scenario files with `dt` that does not divide `t_final`, the output interval, or the controller period, or with events off step boundaries — loader must reject with a clear message (Task 1.4, Task 7.3).
3. Degenerate inputs: zero airspeed (α/β undefined), θ = ±90° Euler extraction, zero rotation-vector exp/log, singular inertia — defined behavior, tested (Tasks 1.1, 4.4).
4. Out-of-table aero lookups during trim iterations or maneuvers — flagged, never silently wrong (Task 4.2).
5. Non-finite state mid-run (NaN from a bad force model) — runner stops and reports time + state (Task 1.4).

---

## Execution topology (orchestrator)

Workers run in separate Paradev workspaces/branches cut from the integration branch; the orchestrator reviews, runs verification, and merges.

| Track | Tasks | Depends on | Suggested agent |
|---|---|---|---|
| A (core) | M0 → M1 → M2 → M3 | — | Claude (Sonnet) for M0/M1/M3; Codex gpt-6.1-sol for M2 math |
| B (spike) | 4.0 JSBSim spike | M0 python env only (can bootstrap its own) | Claude (Sonnet) |
| C (data) | 4.1 F-16 data transcription + reference values | none | Codex gpt-6.1-sol |
| D | M4 model (4.2–4.6) | A through M1, C | Codex |
| E | M5 | D | Codex |
| F | M6 | B, E | Claude/Codex |
| G | M7 | E (and F preferably) | Codex |
| H | M8 | all | Claude |

B and C run in parallel with A from the start. Independent reviews of math-heavy merges by a different provider (Codex reviews Claude work and vice versa; Gemini optional second opinion).

---

## Milestone M0 — Foundations

### Task 0.1: Build system, dependencies, CI, Python env

**Files:**
- Create: `CMakeLists.txt`, `CMakePresets.json`, `cmake/Dependencies.cmake`, `cmake/Warnings.cmake`, `src/flightsim/version.hpp`, `src/flightsim/version.cpp`, `apps/flightsim.cpp`, `tests/unit/test_smoke.cpp`, `tests/validation/test_smoke_validation.cpp`, `bench/bench_smoke.cpp`, `.clang-format`, `.gitignore`, `LICENSE` (MIT, "Eliot Boda"), `README.md` (stub), `.github/workflows/ci.yml`, `python/pyproject.toml`, `python/.python-version` (`3.12`), `python/uv.lock`, `python/flightsim_analysis/__init__.py`, `python/tests/test_smoke.py`

**Interfaces:**
- Produces: target `flightsim_core` (STATIC, PUBLIC include `src/`, links `Eigen3::Eigen`, `tomlplusplus::tomlplusplus`); executable `flightsim`; test executables `unit_tests`, `validation_tests` discovered via `gtest_discover_tests` with `PROPERTIES LABELS unit|validation`; `bench_flightsim`. Presets: `dev` (Debug), `release` (Release), `asan` (Debug, `-fsanitize=address,undefined -fno-omit-frame-pointer`), each with matching build/test presets. `flightsim --version` prints `flightsim <version> (<git sha>)`; git sha injected at configure time.

- [ ] **Step 1:** Write `tests/unit/test_smoke.cpp` asserting `flightsim::version_string()` is non-empty and `Eigen::Vector3d::UnitX().dot(Eigen::Vector3d::UnitX()) == 1.0`.
- [ ] **Step 2:** Write CMake with FetchContent (URL + `URL_HASH SHA256=`; compute hashes with `shasum -a 256`), `FETCHCONTENT_UPDATES_DISCONNECTED ON`, `SYSTEM` on deps, `BENCHMARK_ENABLE_TESTING OFF`, `EIGEN_BUILD_DOC OFF`/`BUILD_TESTING` off for deps. Eigen URL: `https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz`.
- [ ] **Step 3:** Run `cmake --preset release && cmake --build --preset release -j && ctest --preset release` → all pass. Repeat with `dev` and `asan`.
- [ ] **Step 4:** CI `ci.yml`: jobs `build-test` (ubuntu-24.04 × {gcc-13, clang-18} × {dev, release}; macos-14 release), `sanitizers` (ubuntu clang asan preset), `format` (`clang-format --dry-run --Werror` on `src apps tests bench`), `python` (`uv sync --frozen && uv run pytest`), `bench` (release build, run `bench_flightsim --benchmark_min_time=0.05s`, upload JSON artifact). Cache `_deps`.
- [ ] **Step 5:** `cd python && uv init --no-workspace` equivalent → pyproject with deps `numpy scipy plotly kaleido control jsbsim pytest`; `uv lock`; `uv run pytest` passes `test_smoke.py` (imports numpy, jsbsim).
- [ ] **Step 6:** Commit `build: CMake/presets/deps, CI, python env`.

**Accept (M0):** clean clone → `cmake --preset release && cmake --build --preset release && ctest --preset release` passes locally (macOS AppleClang 15); `asan` preset passes; `uv run pytest` passes.

---

## Milestone M1 — Rigid-body core + translational validation

### Task 1.1: Math utilities (units, quaternion, frames)

**Files:** Create `src/flightsim/math/units.hpp`, `src/flightsim/math/quat.hpp`, `src/flightsim/math/quat.cpp`, `src/flightsim/math/frames.hpp`, `src/flightsim/math/frames.cpp`; Test `tests/unit/test_quat.cpp`, `tests/unit/test_frames.cpp`

**Interfaces (Produces):**
```cpp
namespace flightsim {
using Vec3 = Eigen::Vector3d; using Vec4 = Eigen::Vector4d; using Mat3 = Eigen::Matrix3d;
template <int N> using Vec = Eigen::Matrix<double, N, 1>;
namespace units { inline constexpr double ft = 0.3048, slug = 14.59390294, lbf = 4.4482216152605, deg = std::numbers::pi/180.0, g0 = 9.80665; }
Vec4 quat_mul(const Vec4& a, const Vec4& b);          // Hamilton, scalar-first
Vec4 quat_conj(const Vec4& q);
Vec4 quat_normalized(const Vec4& q);
Mat3 dcm_from_quat(const Vec4& q_NB);                 // returns C_NB
Vec3 quat_rotate(const Vec4& q_NB, const Vec3& v_B);  // = C_NB v_B
Vec4 quat_from_euler(double phi, double theta, double psi);   // 3-2-1
Vec3 euler_from_quat(const Vec4& q_NB);               // (phi, theta, psi); theta=±pi/2 → phi=0
Vec4 quat_from_rotvec(const Vec3& rv);                // rotation by |rv| about rv/|rv|; Taylor branch |rv|<1e-4
Vec3 rotvec_from_quat(const Vec4& q);                 // inverse, shortest (w>=0 canonicalization)
double attitude_error(const Vec4& q_ref, const Vec4& q);      // 2*atan2(|vec(q_ref^* ⊗ q)|, |scalar|)
Mat3 skew(const Vec3& v);
struct AirData { double V, alpha, beta; };            // in frames.hpp
AirData air_data_from_body_velocity(const Vec3& v_air_B);    // V<1e-9 → alpha=beta=0
Mat3 euler_rate_matrix(double phi, double theta);     // [phi_dot,theta_dot,psi_dot] = E * omega_B
}
```

- [ ] **Step 1: Failing tests** (`test_quat.cpp`): `HamiltonProductMatchesHandComputed` (i⊗j = k, j⊗i = −k, explicit numeric pair); `ConjugateInverse` (`q⊗q* = [1,0,0,0]` ≤1e-15); `RotateBasis90` (ψ=90°: body x → [0,1,0] (East); θ=90°: body x → [0,0,−1] (Up); φ=90°: body y → [0,0,1]); `DcmOrthonormal` (‖CᵀC−I‖ ≤ 1e-14, det = 1 ± 1e-14 for 1000 seeded random q); `EulerRoundTrip` (grid φ,ψ ∈ [−170°,170°], θ ∈ [−89.9°, 89.9°] step 10°: max error ≤ 1e-12); `GimbalLockDefined` (θ=90° → φ=0 and C reproduced to 1e-12); `ErrorAngleSignInvariant` (`attitude_error(q, −q) == 0`, 1° rotation → 1° ± 1e-12, 1e-9 rad rotation resolved to ≤1e-15 abs); `RotvecRoundTrip` incl. |rv| ∈ {0, 1e-10, 1e-5, 1, 3.1}; `EigenStorageOrderGuard` (`dcm_from_quat(q)` equals `Eigen::Quaterniond(q0,q1,q2,q3).toRotationMatrix()`).
- [ ] **Step 2: Failing tests** (`test_frames.cpp`): `AlphaPositiveForPositiveW`, `BetaPositiveForPositiveV`, `ZeroAirspeedDefined`, `EulerRateMatrixMatchesFiniteDifference` (differentiate `euler_from_quat` along `q̇ = ½ q⊗[0,ω]`, 1e-7).
- [ ] **Step 3:** Run `ctest --preset dev -L unit` → FAIL (undefined). Implement. Run → PASS.
- [ ] **Step 4:** Commit `feat(math): quaternion, frames, units with convention tests`.

### Task 1.2: Rigid-body EOM, forces, integrators (Euler, RK4)

**Files:** Create `src/flightsim/core/state.hpp`, `src/flightsim/core/rigid_body.hpp/.cpp`, `src/flightsim/core/forces.hpp/.cpp`, `src/flightsim/core/integrators.hpp`; Test `tests/unit/test_rigid_body.cpp`, `tests/unit/test_integrators.cpp`

**Interfaces (Produces):**
```cpp
namespace flightsim::rb { inline constexpr int kPos = 0, kVel = 3, kQuat = 6, kOmega = 10, kN = 13; }
// accessors: auto pos(Vec<N>&), vel, quat, omega → Eigen segments (const and non-const), templated on N >= 13
struct MassProps { double mass; Mat3 inertia; Vec3 h_rotor = Vec3::Zero(); Mat3 inertia_inv() const; };
struct ForceMoment { Vec3 F_N = Vec3::Zero(); Vec3 F_B = Vec3::Zero(); Vec3 M_B = Vec3::Zero(); };
ForceMoment& operator+=(ForceMoment&, const ForceMoment&);
// §5.1. g_N = gravity vector in NED. Writes first 13 entries.
Vec<13> rigid_body_deriv(const Vec<13>& x, const MassProps& mp, const ForceMoment& fm, const Vec3& g_N);
using ForceModel = std::function<ForceMoment(double t, const Vec<13>& x)>;
ForceModel linear_drag(double c);           // F_N = -c v_N
ForceModel quadratic_drag(double k);        // F_N = -k |v_N| v_N
ForceModel spring(const Vec3& k_diag);      // F_N = -diag(k) p_N
ForceModel body_force(const Vec3& F_B);
ForceModel body_torque(const Vec3& M_B);
enum class Integrator { Euler, RK4, SymplecticEuler, LieEuler, RKMK4 };
// f: (double t, const Vec<N>& x) -> Vec<N>
template <int N, class F> Vec<N> euler_step(F&& f, double t, const Vec<N>& x, double h);
template <int N, class F> Vec<N> rk4_step(F&& f, double t, const Vec<N>& x, double h);
template <int N> void renormalize_quat(Vec<N>& x);   // N >= 13, layout rb::
struct StepStats { long nfev = 0; };                 // derivative evaluations counter (for work-precision)
```
- [ ] **Step 1: Failing tests** `test_rigid_body.cpp`: `GravityOnlyAcceleration` (v̇ = g_N); `BodyForceRotatedToNed` (q for ψ=90°, F_B = x̂ → v̇_N = [0, F/m, 0]); `QuatKinematicsMatchesFormula` (q̇ = ½ q⊗[0,ω]); `EulerEquationsGyroscopic` (I=diag(1,2,3), ω=(1,1,1) → ω̇ = I⁻¹(−ω×Iω) = (−1, 1, −1/3)); `RotorTermIncluded` (h_rotor=(160,0,0), ω=(0,1,0) → ω̇_z = −(ω×h)_z/Izz).
- [ ] **Step 2: Failing tests** `test_integrators.cpp`: `EulerScalarExponential` (ẋ = −x, one step h=0.1 → 0.9), `Rk4ScalarExponential` (one step error vs e^{−0.1} ≤ 1e-7 and matches the RK4 polynomial `1 − h + h²/2 − h³/6 + h⁴/24` to 1e-15), `RenormalizeLeavesOthers`.
- [ ] **Step 3:** Implement; run unit tests → PASS. Commit `feat(core): rigid-body EOM, force models, Euler/RK4`.

### Task 1.3: Analytic translational references + validation T0–T5

**Files:** Create `src/flightsim/validation/analytic.hpp/.cpp`, `src/flightsim/validation/cases.hpp/.cpp`, `src/flightsim/validation/convergence.hpp/.cpp`; Test `tests/validation/test_translational.cpp`

**Interfaces (Produces):**
```cpp
namespace flightsim::validation {
struct CaseDef {                       // a validation problem with known/fine reference
  std::string name; MassProps mp; Vec3 g_N; std::vector<ForceModel> forces; Vec<13> x0; double t_final;
  std::function<Vec<13>(double t)> exact;   // empty → no closed form
};
CaseDef case_T0(); CaseDef case_T1(); CaseDef case_T2(); CaseDef case_T3(); CaseDef case_T4(); CaseDef case_T5();
struct ErrorNorms { double pos, vel, att, omega; };
ErrorNorms state_error(const Vec<13>& x, const Vec<13>& ref);
struct RunResult { Vec<13> x_final; long nfev; };
RunResult integrate_case(const CaseDef& c, Integrator integ, double h, bool renormalize = true);
struct ConvergenceRow { std::string case_name; Integrator integ; double h; ErrorNorms err; long nfev; };
std::vector<ConvergenceRow> convergence_study(const CaseDef& c, Integrator integ, double h0, int levels);
// least-squares slope over asymptotic window (§6.3): levels with err > 1e3*eps*scale and successive p_k within 10%
struct OrderFit { double order; int first_level, last_level; };
OrderFit observed_order(const std::vector<double>& h, const std::vector<double>& err, double scale);
}
```
Case parameters (fixed, document in code): T2 m=2 kg, c=0.5 N·s/m, v0=(10,−5,−20) m/s, g=9.80665, T=5 s. T3 m=1, k=0.01 kg/m (v_t=√(mg/k)), drop from rest, T=10 s, `z = (v_t²/g) ln cosh(g t/v_t)`. T4 m=1, k=(1,4,9) N/m, p0=(1,1,1), v0=(0,0.5,0), g=0, T=20 s. T5 m=1, I=I₃, ω=(0,0,Ω=2) rad/s constant (torque-free about principal axis), F_B=(1,0,0) N, g=0: `v_N(t) = v0 + (F/(mΩ)) [sin Ωt, 1−cos Ωt, 0]`, p by integration.
- [ ] **Step 1: Failing tests**: `T0ZeroForceExactAllIntegrators` (≤1e-12, smoke), `T1ConstantAccelRk4Exact` (smoke, ≤1e-10), `T2EulerOrder1` / `T2Rk4Order4` (`observed_order` ∈ [0.95,1.05] / [3.9,4.1], h0 = 0.1, levels = 8), `T3Rk4Order4`, `T3EulerOrder1`, `T4EulerEnergyGrows`, `T4Rk4EnergyDriftSmall` (|ΔE|/E ≤ 1e-6 at h=0.01, T=20), `T5DirectionOfRotation` (RK4 h=1e-3 error ≤ 1e-9; and a deliberately transposed-DCM variant computed in-test diverges > 0.1 — documents what the test catches), `ObservedOrderFitExcludesRoundoffFloor` (synthetic data).
- [ ] **Step 2:** Implement; run `ctest --preset release -L validation` → PASS. Commit `test(validation): translational analytic cases T0–T5 and convergence harness`.

### Task 1.4: Scenario loader, runner, CSV/metadata output, CLI `run`

**Files:** Create `src/flightsim/sim/scenario.hpp/.cpp`, `src/flightsim/sim/runner.hpp/.cpp`, `src/flightsim/sim/output.hpp/.cpp`, `scenarios/rb_linear_drag.toml`, `scenarios/rb_torque_free_asym.toml`, `scenarios/rb_dzhanibekov.toml`; Modify `apps/flightsim.cpp`; Test `tests/unit/test_scenario.cpp`, `tests/unit/test_runner.cpp`

**Interfaces (Produces):**
```cpp
namespace flightsim::sim {
struct RunConfig { Integrator integ; double dt, t_final; int output_every = 1; bool renormalize = true; };
struct RigidBodyScenario { RunConfig run; MassProps mp; Vec3 g_N; std::vector<ForceModel> forces; Vec<13> x0; std::string source_path, source_hash; };
RigidBodyScenario load_rigid_body_scenario(const std::filesystem::path&);   // throws std::runtime_error with key path in message
enum class StopReason { FinalTime, GroundContact, NonFinite };
struct Sample { double t; Vec<13> x; };
struct RunOutput { std::vector<Sample> samples; StopReason reason; double t_stop; long nfev; };
RunOutput run(const RigidBodyScenario&);
void write_csv(const std::filesystem::path&, const RigidBodyScenario&, const RunOutput&);   // columns per §9.3 incl. derived (phi_rad..., E_J, HN_x..., qnorm)
void write_metadata(const std::filesystem::path&, const RigidBodyScenario&, const RunOutput&); // JSON via toml::json_formatter
}
```
TOML schema = PRD §9.4 (`[sim]`, `[body]`, `[initial]`, `[[forces]] type = "linear_drag"|"quadratic_drag"|"spring"|"body_force"|"body_torque"` + params, `[environment] gravity_mps2 = 9.80665`). Unknown keys → error.
- [ ] **Step 1: Failing tests**: `RejectsDtNotDividingTFinal`, `RejectsUnknownKey`, `RejectsSingularInertia`, `RejectsMissingField`, `LoadsExampleScenario`, `NonFiniteStops` (force returning NaN at t≥0.5 → reason NonFinite, t_stop = first step boundary ≥ 0.5, message lists state), `GroundContactStopsWhenHBelowZero` (only when `[sim] stop_at_ground = true`), `DeterministicBitwise` (two runs → identical CSV bytes).
- [ ] **Step 2:** Implement + CLI `flightsim run <scenario> --out <dir>` (writes `trajectory.csv`, `metadata.json`). Run tests → PASS. Run `./build/release/apps/flightsim run scenarios/rb_linear_drag.toml --out out/rb_linear_drag` → files exist.
- [ ] **Step 3:** Commit `feat(sim): TOML scenarios, runner, CSV/JSON output, CLI run`.

**Accept (M1):** T2/T3 orders in bounds; T5 passes; every §4 convention has a test (checklist in PR description); determinism test passes; `asan` preset green.

---

## Milestone M2 — Rotational validation, convergence harness, extra integrators (MVP)

### Task 2.1: Symplectic Euler, Lie–Euler, RKMK4

**Files:** Modify `src/flightsim/core/integrators.hpp`; Test `tests/unit/test_lie_integrators.cpp`

**Interfaces (Produces):** `symplectic_euler_step(f, t, x, h)` (N ≥ 13: `v⁺ = v + h·f(x)_vel; p⁺ = p + h·v⁺`; quat/ω use Euler — documented non-symplectic), `lie_euler_step`, `rkmk4_step` (template N ≥ 13, rb layout; attitude via §6.1 formula **Θ̇ = ω + ½Θ×ω + (1/12)Θ×(Θ×ω)**, `q⁺ = q ⊗ quat_from_rotvec(Θ)`, other states RK4; quaternion component of `f` unused). `step(Integrator, f, t, x, h, bool renormalize, StepStats&)` dispatcher (renormalize applies to Euler/RK4/SymplecticEuler only).
- [ ] **Step 1: Failing tests**: `LieEulerPreservesNorm` (10⁶ steps, |‖q‖−1| ≤ 1e-13), `Rkmk4ExactForConstantOmega` (R1: ≤1e-12 rad after 1000 steps), `Rkmk4SignRegression` (one step on a case with Θ ∦ ω — compare with tiny-step RK4 reference: error ≤ C·h⁵ with h=0.01 → ≤ 1e-9; the "−½" variant computed in-test is ≥ 100× worse), `SymplecticEulerOscillatorEnergyBounded` (T4, 10⁵ steps h=0.01: max |ΔE|/E ≤ 0.05 and no secular trend: |mean over last 10% − mean over first 10%| ≤ 1e-3).
- [ ] **Step 2:** Implement; PASS; commit `feat(core): symplectic Euler, Lie–Euler, RKMK4 integrators`.

### Task 2.2: Rotational references R1–R5, Jacobi elliptic functions

**Files:** Create `src/flightsim/math/elliptic.hpp/.cpp`, `tests/data/ellipj_reference.csv`, `python/tools/gen_ellipj_reference.py`; Modify `src/flightsim/validation/analytic.*`, `cases.*`; Test `tests/unit/test_elliptic.cpp`, `tests/validation/test_rotational.cpp`

**Interfaces (Produces):** `struct JacobiSnCnDn { double sn, cn, dn; }; JacobiSnCnDn jacobi_elliptic(double u, double m);` (parameter m = k², AGM A&S 16.4, valid 0 ≤ m < 1). `CaseDef case_R1(), case_R2(), case_R3(), case_R4(), case_R5();` Conservation helpers: `double rotational_energy(const Vec<13>&, const MassProps&)`, `Vec3 angular_momentum_N(const Vec<13>&, const MassProps&)` (incl. rotor).
Parameters: R1 I=diag(1,2,3), ω=(0,0,1.3), T=10. R2 I=diag(2,2,1), ω0=(0.3,0.4,2.0), q0 from Euler(10°,20°,30°), T=10, formula §8.3 with ĥ_N = H_N/|H_N|. R3 I=diag(1,2,3), ω0=(0.5,0,1.0) (L²=9.25 > 2E·I₂=6.5), Landau–Lifshitz §37 formulas with τ = t·sqrt((I3−I2)(L²−2E I1)/(I1 I2 I3)), m = (I2−I1)(2E I3−L²)/((I3−I2)(L²−2E I1)), T=20. R4 I=diag(1,2,3), ω0=(1e-3, 1.0, 1e-3), T=60 (demo). R5 I=diag(1,2,3), h_rotor=(0.5,0,0), ω0=(0.2,0.3,0.4), T=20.
- [ ] **Step 1:** `gen_ellipj_reference.py` (scipy.special.ellipj, grid u ∈ [−5,5], m ∈ {0, 0.1, 0.5, 0.9, 0.99, 0.999}); commit CSV.
- [ ] **Step 2: Failing tests** `test_elliptic.cpp`: `MatchesScipy` (≤1e-12 abs), `Identities` (sn²+cn²=1, dn²+m sn²=1 ≤1e-14), `LimitsM0AndM1` (m=0 → sin/cos/1).
- [ ] **Step 3: Failing tests** `test_rotational.cpp`: `R2ReferenceSatisfiesKinematics` (central FD of formula vs ½q⊗ω ≤ 1e-9), `R3ReferenceSatisfiesEuler` (FD residual ≤ 1e-9; conserves E, |L| ≤ 1e-13), `R1Rk4Order4`, `R2Rk4Order4`, `R2Rkmk4Order4`, `R2LieEulerOrder1`, `R2EulerOrder1`, `R3*` same set on ω-error and on H_N drift, `HNConservedRk4` (R3 h=1e-3: ‖ΔH_N‖/‖H_N‖ ≤ 1e-10), `R5RotorHNConserved`, `NormDriftWithoutRenormalization` (Euler/RK4 renormalize=false reported; RKMK4 ≤ 1e-13), `R4FlipOccursAndConserves` (ω₂ changes sign at least once within T; E, |H| ≤ 1e-8 at h=1e-3 RK4).
- [ ] **Step 4:** Implement, PASS, commit `test(validation): rotational references R1–R5 with elliptic functions`.

### Task 2.3: Generic linearizer + L1; CLI `converge`

**Files:** Create `src/flightsim/analysis/jacobian.hpp`; Modify `apps/flightsim.cpp`, `src/flightsim/validation/convergence.*`; Create `scenarios/studies/core_convergence.toml`; Test `tests/unit/test_jacobian.cpp`, `tests/validation/test_spin_stability.cpp`

**Interfaces (Produces):** `template <int N, int M, class F> Eigen::Matrix<double,N,M> jacobian_central(F&& f, const Vec<M>& x, const Vec<M>& scale, double eps = 6e-6);` (δ_i = eps·max(|x_i|, scale_i)). CLI `flightsim converge <study.toml> --out <dir>` → `convergence.csv` (case, integrator, h, err_pos, err_vel, err_att, err_omega, err_E, err_H, nfev) + `orders.json`. Study TOML: `cases = ["T2",...]`, `integrators = [...]`, `h0`, `levels` per case.
- [ ] **Step 1: Failing tests**: `JacobianOfLinearMapExact` (≤1e-9), `L1SpinStability` (Euler's equations ω̇ = I⁻¹(−ω×Iω), linearized at ω = Ω e_i for each axis i, I=diag(1,2,3), Ω=1.5: eigenvalues vs λ² = −Ω²(Ii−Ij)(Ii−Ik)/(Ij Ik) ≤ 1e-6 rel; intermediate axis has a positive real eigenvalue).
- [ ] **Step 2:** Implement, PASS. Run `flightsim converge scenarios/studies/core_convergence.toml --out out/convergence` → orders.json matches PRD §10.3 table. Commit.

**Accept (M2 = MVP):** all §10.3 order assertions pass for all five integrators; RK4 vs RKMK4 norm/conservation comparison emitted; L1 passes; full suite green on dev/release/asan.

---

## Milestone M3 — Analysis tooling, demo, docs, benchmarks (resume-ready core)

### Task 3.1: Benchmarks + zero-allocation test
**Files:** `bench/bench_integrators.cpp`, `tests/alloc/test_no_alloc.cpp` (own executable `alloc_tests`, label `unit`, global `operator new` counter). Benchmarks: rigid-body RK4 step (gravity + quadratic drag), Euler/RK4/RKMK4 step cost, R3 10 s run per integrator. Accept: alloc count = 0 over 10⁴ RK4 steps after warmup.

### Task 3.2: Python analysis + report (use **Impeccable** for the HTML report design)
**Files:** `python/flightsim_analysis/{run.py,plots.py,report.py,metrics.py}`, `python/tests/test_report.py`, `scripts/demo.sh`, `docs/results.template.md`.
- `run.py`: invokes CLI on every `scenarios/*.toml` and studies; collects `out/`.
- `plots.py`: Plotly figures — convergence log-log with reference slopes, work–precision, energy/H drift vs time per integrator, quaternion norm drift, R4 Dzhanibekov ω(t) + H_N components, T4 phase portrait. Labeled axes with units; HTML + PNG (kaleido).
- `report.py`: writes `out/index.html` (static, single page, figures embedded) and regenerates `docs/results.md` from `docs/results.template.md` + `out/metrics.json` (no hand-typed numbers).
- `scripts/demo.sh`: `cmake --preset release && cmake --build --preset release -j && ctest --preset release && (cd python && uv sync --frozen && uv run python -m flightsim_analysis.run && uv run python -m flightsim_analysis.report)`; prints path to `out/index.html`.
- Test: `test_report.py` regenerates results.md from a fixture metrics.json and asserts every `{{key}}` is filled.

### Task 3.3: Docs
**Files:** `docs/design.md` (architecture, conventions §4 restated with test names, EOM/quaternion/RKMK4 derivations, decision log, limitations), `README.md` (pitch, headline plots, numbers from results.md, quickstart). Accept (M3): demo from clean clone < 5 min excluding dependency download; regenerates every core result.

---

## Milestone M4 — F-16 model

### Task 4.0: JSBSim feasibility spike (throwaway, timeboxed; parallel track B)
**Files:** `spikes/jsbsim/` (throwaway scripts, not product code), findings in `docs/design.md#jsbsim-spike` (create section) and `PROGRESS.md`.
Answer with evidence (script + output) each of: (1) `jsbsim.FGFDMExec(root_dir)` loads a custom aircraft from a custom directory; (2) a custom planet with `rotation_rate = 0` loads from Python (`fdm.load_planet(...)` or equivalent; record exact API / version); (3) `simulation/gravity-model` settable to standard (0); (4) integrator properties `simulation/integrator/rate/rotational|translational`, `simulation/integrator/position/rotational|translational` settable, list allowed values; (5) a minimal custom aircraft with body-axis aero functions and an external force flies and can be trimmed; (6) surfaces can be driven directly (no FCS) in a custom aircraft; (7) stock `f16`: can surface positions be set directly bypassing its FCS? Its data provenance (read `aircraft/f16/f16.xml` header). (8) `<mass_balance negated_crossproduct_inertia="false">` accepted and reflected in `inertia/ixz-slugs_ft2`. Output: a short table yes/no + exact API, and a recommendation on whether PRD §8.6 needs the fallback.

### Task 4.1: F-16 data transcription + reference values (parallel track C)
**Files:** `data/f16/f16_sl.toml`, `python/tools/f16_data_diff.py`, `data/f16/SOURCES.md`, `data/f16/reference_values.toml`
- Transcribe every S&L low-fidelity table and constant (Cx(α,δe), Cz(α) + β/δe terms, Cm(α,δe), Cl(α,β), Cn(α,β), DLDA, DLDR, DNDA, DNDR, damping table (9 coefficients vs α), thrust idle/mil/max (alt × Mach), throttle gearing, power-lag rtau, mass/geometry, x_cg_ref, h_x) **numbers only** from source A (AeroBenchVV, MATLAB or Python), with breakpoints and units.
- Independently transcribe from source B (F16Model.jl or Archimedes) into a scratch file; `f16_data_diff.py` diffs A vs B value-by-value; resolve disagreements against TP-1538 (NTRS 19800005879) where possible; record every resolution in `SOURCES.md`.
- `reference_values.toml`: A1 check case inputs and expected derivatives, TR1 trim results, LN1 linear model / eigenvalues — each value with ≥ 2 sources or flagged `single_source = true`; note disagreements.
- Accept: diff script exits 0 (or only documented resolutions); every table has a citation; licenses of A and B recorded.

### Task 4.2 (detailed): Interpolation + F-16 data loader

**Files:** Create `src/flightsim/math/interp.hpp/.cpp`, `src/flightsim/aircraft/f16_data.hpp/.cpp`; Test `tests/unit/test_interp.cpp`, `tests/unit/test_f16_data.cpp`.

**Interfaces (Produces):**
```cpp
namespace flightsim {
struct LookupMonitor {                 // diagnostics; not thread-safe; passed by pointer, may be null
  bool out_of_domain = false; std::string first_violation;   // "<table>: <axis>=<value> outside [lo,hi]"
  bool record_cells = false; std::vector<std::array<int,2>> cells;   // per lookup (table-id hash, cell index) when record_cells
  void reset();
};
class Table1D { public: Table1D(std::string name, std::vector<double> x, std::vector<double> y);
  double operator()(double x, LookupMonitor* mon = nullptr) const; };
class Table2D { public: Table2D(std::string name, std::vector<double> x, std::vector<double> y, Eigen::MatrixXd v);  // v(i,j) at (x_i, y_j)
  double operator()(double x, double y, LookupMonitor* mon = nullptr) const; };
// Semantics (S&L): piecewise (bi)linear; outside the breakpoint range, linear extrapolation using the boundary cell
// (no distance limit), and the monitor flags out_of_domain. Breakpoints strictly increasing (checked in ctor).
struct F16Data { /* all SI; angles rad; see header doc for every field */ };
F16Data load_f16_data(const std::filesystem::path& toml_path);   // throws std::runtime_error naming the bad key
}
```
F16Data rules: aero coefficient tables stay dimensionless; angle breakpoints converted deg→rad with `units::deg`; mass = `1/reciprocal_mass` slug → kg (S&L uses RM, not W/g — document); inertia slug·ft² → kg·m² (Ixz stored as the positive product ∫xz dm per S&L); S, b, cbar → m², m; `hx` → kg·m²/s; thrust tables lbf → N with altitude breakpoints ft → m, Mach dimensionless; surface limits deg → rad; actuator lag 1/s and rate limits → rad/s; engine throttle gearing / rtau / transition constants unchanged (percent, 1/s). Every TOML table maps to one named field; unknown or missing tables are errors. The header documents each field's meaning, units, axis order and source key.
- [ ] Tests: `ExactAtBreakpoints`, `LinearBetweenBreakpoints`, `ExtrapolatesFromBoundaryCell` (both ends, 1D and 2D), `MonitorFlagsOutOfDomain` (+ message names table and axis), `MonitorRecordsCells`, `RejectsNonIncreasingBreakpoints`, `RejectsShapeMismatch`; `LoadsRealDataFile` (counts: cx 12×5, cl 12×7, damping 9 tables, thrust 6×6), `MassPropsConvertedToSI` (mass = (1/0.00157)·14.59390294 kg ± 1e-9 rel; Ixx = 9496·14.59390294·0.3048² ± 1e-12 rel), `AngleBreakpointsInRadians`, `ThrustInNewtons`, `RejectsMissingTable` (temp copy with one table removed).

### Task 4.3 (detailed): Atmosphere + gravity

**Files:** Create `src/flightsim/env/atmosphere.hpp/.cpp`; Test `tests/unit/test_atmosphere.cpp`.

**Interfaces (Produces):**
```cpp
namespace flightsim {
enum class AtmosphereKind { SlAdc, Us76 };
struct AtmState { double rho, T, p, a; };            // kg/m³, K, Pa, m/s
AtmState atmosphere(AtmosphereKind, double h_m);      // geometric altitude; valid 0..20 km (Us76), S&L range for SlAdc; throws std::domain_error outside [-500 m, 20 km]
struct GravityModel { double g0 = units::g0; double R = std::numeric_limits<double>::infinity();
  double gravity(double h_m) const; };              // g0·(R/(R+h))²; R=inf → g0
}
```
- SlAdc = S&L `adc` routine: `tfac = 1 − 0.703e-5·h_ft`, `T = 519·tfac °R` (390 °R for h ≥ 35,000 ft), `rho = 2.377e-3·tfac^4.14 slug/ft³`, `a = sqrt(1.4·1716.3·T) ft/s`, converted to SI. Constants hard-coded in `atmosphere.cpp` with a test (`SlAdcConstantsMatchDataFile`) that they equal `data/f16/f16_sl.toml [environment]`. p from the ideal gas law with the same gas constant.
- Us76 = U.S. Standard Atmosphere 1976 layers 0–11 km (−6.5 K/km) and 11–20 km (isothermal), in **geopotential** altitude `H = r0·h/(r0+h)`, r0 = 6,356,766 m; T0 = 288.15 K, p0 = 101,325 Pa, R = 287.05287 J/(kg·K), g0 = 9.80665, γ = 1.4.
- [ ] Tests: `Us76MatchesStandardTable` at 0, 1, 5, 11, 15, 20 km for ρ, T, p, a ≤ 1e-4 rel against values transcribed (with citation in the test) from the published 1976 standard tables; `SlAdcMatchesRoutine` (hand-evaluated formula at 0, 10,000, 40,000 ft); `SlAdcConstantsMatchDataFile`; `GravityInverseSquare` (R = 6,371,000 m, h = 10 km); `GravityConstantWhenRInfinite`; `ThrowsOutsideRange`.

### Task 4.4 (detailed): Wind model

**Files:** Create `src/flightsim/env/wind.hpp/.cpp`; Test `tests/unit/test_wind.cpp`.
**Interfaces (Produces):**
```cpp
namespace flightsim {
struct Gust { double t_start_s; double duration_s; double amplitude_mps; Vec3 direction_N; };   // full-wave 1-cos
struct WindModel { Vec3 constant_N = Vec3::Zero(); std::optional<Gust> gust;
  Vec3 wind_N(double t) const; };   // constant + gust: A/2·(1 − cos(2π (t−t0)/T)) for t0 ≤ t ≤ t0+T, else 0; direction normalized
}
```
Shape choice (document in header + design.md): MIL-F-8785C discrete-gust 1-cos profile, applied rise-and-fall over `duration` (frozen in time, spatial gradients neglected per PRD §5.4).
- [ ] Tests: `ConstantOnly`, `GustZeroBeforeAndAfter`, `GustPeakAtMidpoint` (= amplitude·dir ± 1e-15), `GustContinuous` (value at t0 and t0+T is 0), `DirectionNormalized`.

### Task 4.5 (detailed): F-16 model + flight-variable chart

**Files:** Create `src/flightsim/aircraft/f16_model.hpp/.cpp`, `src/flightsim/analysis/flight_vars.hpp/.cpp`; Test `tests/unit/test_f16_model.cpp`, `tests/unit/test_flight_vars.cpp`.
**Interfaces (Produces):**
```cpp
namespace flightsim {
namespace ac { inline constexpr int kPower = 13, kDe = 14, kDa = 15, kDr = 16, kN = 17; }
struct F16Input { double throttle = 0, de = 0, da = 0, dr = 0; };   // throttle [0,1]; surfaces rad (S&L signs)
struct F16Model {
  F16Data data; double xcg = 0.35;
  AtmosphereKind atm = AtmosphereKind::SlAdc; GravityModel grav;   // caller sets grav.g0 (S&L: data.g_sl)
  WindModel wind; bool actuators = false; bool include_rotor = true;
  LookupMonitor* monitor = nullptr;                                 // diagnostics only
  struct Eval { ForceMoment fm; double V, alpha, beta, mach, qbar, thrust, power_cmd; Vec<4> surfaces; };
  Eval evaluate(double t, const Vec<17>& x, const F16Input& u) const;   // aero + thrust, non-gravitational
  Vec<17> deriv(double t, const Vec<17>& x, const F16Input& u) const;   // rigid_body_deriv on first 13 + engine + actuators
  MassProps mass_props() const;     // inertia tensor FRD with off-diagonals −Ixz; h_rotor = (hx,0,0) if include_rotor
};
// x13 = [V, alpha, beta, phi, theta, psi, p, q, r, pN, pE, h, P]  (S&L order; SI; V air-relative)
Vec<17> x17_from_x13(const F16Model&, const Vec<13>& x13, const F16Input& u, double t = 0);   // actuators set to u
Vec<13> x13_from_x17(const F16Model&, const Vec<17>& x, double t = 0);
Vec<13> f13(const F16Model&, const Vec<13>& x13, const F16Input& u, double t = 0);   // assumes steady wind
}
```
Model content = S&L Appendix A as described in `data/f16/SOURCES.md` and the `F16Data` header: Cx, Cy, Cz, Cl+DLDA·δa/norm+DLDR·δr/norm, Cm, Cn+…; damping added with `q·cbar/(2V)` and `p,r·b/(2V)`; CG transfer `Cm += Cz_total·(xcg_ref − xcg)`, `Cn −= Cy_total·(xcg_ref − xcg)·cbar/b` (totals after damping); forces `q̄S·C`, moments `q̄Sb·Cl`, `q̄S·cbar·Cm`, `q̄Sb·Cn`; thrust along +x_B from the power-level interpolation between idle/mil/max tables at (h, Mach); engine `Ṗ` with throttle gearing, rtau and the 50 % transition logic; actuators when enabled: `δ̇ = clamp(lag·(clamp(cmd, limits) − δ), ±rate)`, else surfaces = commands and actuator-state derivatives = 0. Air data from `v_air,B = C_BN (v_N − wind_N(t))`. Gravity `g_N = (0,0,grav.gravity(h))`. `f13` maps body accelerations to `V̇, α̇, β̇` (S&L formulas) and body rates to Euler rates via `euler_rate_matrix`.
- [ ] Tests (`test_flight_vars.cpp`): `RoundTripX13`, `F13EulerRatesMatchFiniteDifference` (vs FD of `x13_from_x17` along `deriv`, 1e-7).
- [ ] Tests (`test_f16_model.cpp`): `A1CheckCase` — model with xcg 0.4, SlAdc, `grav.g0 = data.g_sl`, R = ∞, rotor on, actuators off; inputs from `data/f16/reference_values.toml [A1]` (load in test with toml++); compare all 13 `f13` components (converted back to ft, ft/s, percent) with `[A1.expected]` at ≤ 1e-3 rel, print each relative error, and mark any > 1e-4 with `ADD_FAILURE() << "investigate"` disabled behind a documented list (the Archimedes reproduction itself differs by ~8e-5 on ẇ); also compare `[A1.body_axes]` u̇,v̇,ẇ. Sign tests: `PositiveElevatorNoseDown` (ΔCm < 0), `PositiveAileronLeftRoll` (ΔL < 0), `PositiveRudderLeftYaw` (ΔN < 0), `PitchDampingNegative` (∂Cm/∂q < 0 at α = 5°), `ThrustIncreasesWithThrottle`, `EnginePowerLagTransitions` (P approaches command; crossing 50 % uses targets), `ActuatorRateAndPositionLimits`, `ZeroAirspeedDoesNotNaN` (V = 0 → finite derivative, documented behavior), `MonitorFlagsOutOfTableAlpha` (α = 50°).

### Task 4.6 (detailed): Aircraft consistency laws + convergence

**Files:** Test `tests/validation/test_f16_consistency.cpp`.
Common run: x0 at h = 3048 m (10,000 ft), V = 152.4 m/s (500 ft/s), α = 0.0594 rad, θ = α, throttle 0.1678, δe = −0.6531° (near-trim values from the Archimedes tutorial; not a reference), xcg 0.35, SlAdc, `g0 = data.g_sl`; elevator doublet ±1° (1–2 s, 2–3 s) and aileron doublet ±2° (5–6 s, 6–7 s), events on step boundaries; RK4. Every run asserts `monitor.out_of_domain == false`.
- [ ] `B1WorkEnergy`: augmented state `[x17 | W_nc]`, `Ẇ_nc = (C_NB F_B)·v_N + M_B·ω_B`, `E = ½m|v_N|² + ½ωᵀIω + m·g0·h`: 60 s at dt = 0.01 → `|ΔE − W_nc|/max|ΔE|` < 1e-6; also at dt = 0.02 / 0.04 shows decrease with observed order printed.
- [ ] `B2AngularImpulse`: augmented `[x17 | J_N]`, `J̇_N = C_NB M_B`; `H_N = C_NB(Iω + h_rotor)`; `‖ΔH_N − J_N‖/max‖ΔH_N‖` < 1e-6 at dt = 0.01.
- [ ] `S1MirrorSymmetrized`: test-only copy of data with DLDA/DLDR/DNDA/DNDR replaced by their β-even part `(T(α,β)+T(α,−β))/2`, `include_rotor = false`, ±2° aileron doublets from wings-level: lateral states (β, p, r, φ, ψ, pE) exact negatives and longitudinal identical ≤ 1e-12 rel over 10 s. `S1ProductionAsymmetryReported`: same with real data and rotor on → print max asymmetry (no assertion).
- [ ] `W1GalileanInvariance`: constant wind (10, −5, 0) m/s; v_N0 shifted by the wind; air-relative states (V, α, β, p, q, r, φ, θ, ψ, h) equal the no-wind run ≤ 1e-9 rel over 20 s; ground track offset = w·t ≤ 1e-6 m.
- [ ] `AircraftConvergenceMeasured`: RK4 with h = 0.04 … 0.0025 against an h/16 RK4 reference over 10 s of the doublet run; print observed order on position and attitude; assert 0.9 ≤ order ≤ 4.2 only (PRD §6.3: reduced order expected from non-smooth tables; explanation goes in results.md).

**Accept (M4):** PRD §13 M4 — A1 within tolerance; B1/B2/S1/W1 pass; aircraft convergence order measured and printed.

### (superseded heading) Tasks 4.4–4.6 interface sketch
- 4.2 `aircraft/f16_data.{hpp,cpp}`: `F16Data load_f16_data(path)` (SI-converted), `math/interp.hpp`: `Table1D`, `Table2D` with S&L one-cell linear extrapolation, `LookupMonitor` (out-of-domain flag + optional breakpoint-cell record). Tests: exact at breakpoints, linear between, extrapolation, monitor.
- 4.3 `env/atmosphere.{hpp,cpp}`: `enum class AtmosphereKind { SlAdc, Us76 }`, `struct AtmState { double rho, T, p, a; }`, `AtmState atmosphere(AtmosphereKind, double h_m)`; `GravityModel { double g0; double R = inf; double gravity(double h) const; }`. Tests per §10.2.
- 4.4 `aircraft/f16_model.{hpp,cpp}`: `struct F16Input { double throttle, de, da, dr; }` (rad), `struct F16Model { F16Data data; AtmosphereKind atm; GravityModel grav; WindModel wind; bool actuators; double xcg; Vec<17> deriv(double t, const Vec<17>& x, const F16Input& u) const; ForceMoment forces(...) const; }`; `analysis/flight_vars.{hpp,cpp}`: x13 ↔ model-state maps, `Vec<13> f13(const F16Model&, const Vec<13>& x13, const F16Input& u)`.
- 4.5 Validation: A1, sign tests, S1 (symmetrized copy), W1, B1, B2 (augmented-state budget), aircraft convergence measurement.
- 4.6 Wind: `env/wind.{hpp,cpp}` constant + 1-cos gust (onset aligned to dt, validated at load).
**Accept (M4):** PRD §13 M4.

## Milestone M5 — Trim and linearization

### Task 5.1 (detailed): Aircraft scenario + trim solver
**Files:** Create `src/flightsim/analysis/trim.hpp/.cpp`, `src/flightsim/sim/aircraft_scenario.hpp/.cpp`, `scenarios/f16_trim_502_xcg035.toml`, `scenarios/f16_trim_502_xcg030.toml`; Modify `apps/flightsim.cpp`; Test `tests/unit/test_trim.cpp`, `tests/validation/test_f16_trim.cpp`.
- Aircraft scenario TOML: `[aircraft] data = "data/f16/f16_sl.toml"` (path relative to the repo root or the scenario file), `xcg`, `atmosphere = "sl_adc"|"us76"`, `gravity = "sl"|"standard"|{g0, R}`, `rotor = true`, `actuators = false`; `[trim] V_mps | V_fps, h_m | h_ft, gamma_deg`. Unknown keys are errors. `F16Model load_aircraft(const toml::table&, base_dir)` plus `TrimSpec`.
- Trim: unknowns `z = [α, throttle, δe]`, with β = φ = p = q = r = δa = δr = 0, θ = α + γ, ψ = 0 and P = power_cmd(throttle) (steady engine). Residual `r = [V̇/g, α̇, q̇]` from `f13`. Levenberg–Marquardt with a central-difference Jacobian (`jacobian_central`), λ adapted ×10 / ÷10, at most 50 iterations, converged when ‖r‖∞ < 1e-10. Initial guess α = 5°, throttle 0.2, δe = −1°. Fails (ok = false, message) if throttle leaves [0,1], δe leaves its limits, α leaves the table domain, or there is no convergence.
- CLI `flightsim trim <scenario> --out trim.json` (α, throttle, δe in SI and in S&L units, residual, iterations, the full x17 state).
- [ ] Tests: `TrimConvergesOnSyntheticQuadratic` (solver on a 3-D test function); `TR1Level35`, `TR1Level30`, `TR1Level38` against `reference_values.toml` (α, throttle, δe within ≤ 1 in the last printed digit or ≤ 0.5 % rel, whichever is larger; print the differences); `TR2Stationarity` (10 s RK4 dt = 0.01 from TR1 level_35: |ΔV| < 1e-3 ft/s, |Δα|, |Δθ| < 1e-6 rad, |Δh| < 1e-2 ft; if an unstable mode makes this impossible from a 1e-10 residual, derive the bound from the measured eigenvalue and document it); `TrimReportsFailureOutsideEnvelope` (V = 60 ft/s).

### Task 5.2 (detailed): Linearization + modes
**Files:** Create `src/flightsim/analysis/linearize.hpp/.cpp`, `src/flightsim/analysis/modes.hpp/.cpp`; Modify `apps/flightsim.cpp`; Test `tests/validation/test_f16_linear.cpp`.
- `LinearModel linearize(const F16Model&, const TrimResult&, double eps = 6e-6)`: A = ∂f13/∂x13 (13×13), B = ∂f13/∂u (13×4, u = [throttle, δe, δa, δr] in rad), via `jacobian_central` with per-variable scales (V: 1 m/s, angles: 1 rad, rates: 1 rad/s, positions: 1 m, P: 1 %).
- `struct Mode { std::string name; std::complex<double> eig; double wn, zeta, tau; }`; `std::vector<Mode> identify_modes(const LinearModel&)`: longitudinal subset [V, α, θ, q], lateral subset [β, φ, p, r]; short period = longitudinal complex pair with the larger ωn, phugoid = the smaller; Dutch roll = lateral complex pair; roll = fastest lateral real; spiral = slowest lateral real.
- CLI `flightsim linearize <scenario> --out linear.json` (A, B, state/input names and units, modes).
- [ ] Tests: `LN1Level30Eigenvalues`: xcg = 0.30 scenario vs `[LN1.level_30]` textbook phugoid, short period, Dutch roll, spiral (and roll if present) within 1 % or 1 in the last printed digit. Note the reference's longitudinal values come from a model without the power state ("tabulated engine"): compare the 4-state longitudinal subset with the engine power held at trim, and say so in the test. `LN1LateralMatrix` vs `[LN1.level_30.lateral_textbook].A` (S&L units; ≤ 1 % or 1 in the last digit per entry). `LN1LongitudinalTutorialMatrix` vs `[longitudinal_tutorial]` (reported, labelled "tutorial-computed", asserted at 2 % since it is rounded to display precision). `LN2QuadraticScaling` (Task 4.x monitor; trailing-edge-up δ = 0.5° vs 0.05°, 5 s; asserts no breakpoint cell change; ratio in [70, 130]). `LN3StepSweep` (ε = 1e-2 … 1e-10; prints the max entry change; asserts the default ε lies in the flat region: relative change ≤ 1e-5 versus ε×10 and ε/10).

### Task 5.3 (detailed): Report extension
**Files:** Modify `python/flightsim_analysis/{run.py,metrics.py,plots.py}`, `docs/results.template.md`, `python/tests/fixtures/metrics.json`.
- Run `flightsim trim`/`linearize` on the scenarios; results.md gains "Trim (TR1)" and "Linear model and modes (LN1)" tables with ours vs published values and relative differences, each labelled **single-lineage secondary source** per PRD §8.5; LN3 sweep figure; LN2 linear-vs-nonlinear overlay figure.

**Accept (M5):** PRD §13 M5.

## Milestone M5 — (original sketch, superseded by Tasks 5.1–5.3) Interfaces: `analysis/trim.{hpp,cpp}`: `struct TrimSpec { double V, h, gamma, xcg; }; struct TrimResult { bool ok; double alpha, throttle, de; Vec<17> x; F16Input u; int iterations; double residual; std::string message; }; TrimResult trim_level(const F16Model&, const TrimSpec&);` `analysis/linearize.{hpp,cpp}`: `struct LinearModel { Eigen::Matrix<double,13,13> A; Eigen::Matrix<double,13,4> B; Vec<13> x0; Vec<4> u0; }; LinearModel linearize(const F16Model&, const TrimResult&, double eps = 6e-6);` `analysis/modes.{hpp,cpp}` mode identification. CLI `trim`, `linearize` (JSON). Tests TR1, TR2, LN1, LN2 (breakpoint monitor asserts no crossing), LN3. **Accept:** PRD §13 M5.

## Milestone M6 — JSBSim matched cross-validation

### Task 6.1: Open-loop aircraft runner + CLI
**Files:** Create `src/flightsim/sim/aircraft_runner.hpp/.cpp`; Modify `src/flightsim/sim/aircraft_scenario.*`, `src/flightsim/sim/output.*`, `apps/flightsim.cpp`; Test `tests/unit/test_aircraft_runner.cpp`.
- `flightsim run` dispatches on the presence of an `[aircraft]` table. The aircraft scenario adds: `[initial] from_trim = true` (uses `[trim]`) or explicit x13; `[[inputs]] channel = "throttle"|"de"|"da"|"dr", t_start_s, t_end_s, value (deg for surfaces), mode = "add"` (added to trim, held constant within steps; event times must be multiples of dt → load error otherwise); `[wind]` constant/gust (Task 4.4); `[sim] integrator, dt_s, t_final_s, output_every, stop_at_ground = true`.
- `struct AircraftRunOutput { std::vector<double> t; std::vector<Vec<17>> x; std::vector<F16Input> u; StopReason reason; bool out_of_domain; std::string first_violation; }`; CSV columns: t_s, x17 (SI), x13-derived (V_mps, alpha_rad, beta_rad, phi/theta/psi_rad, h_m), mach, qbar_Pa, thrust_N, throttle, de_rad, da_rad, dr_rad, plus the same in S&L units with `_ft`, `_fps`, `_deg` suffixes for comparison convenience.
- `struct ControllerHook` is NOT added here (M7 adds controller scheduling); keep the loop structure so M7 inserts a fixed-rate controller call before each step.
- [ ] Tests: `RejectsUnalignedInputEvent`, `DoubletAppliedOnSchedule`, `FromTrimStaysTrimmedWithNoInputs` (2 s, |ΔV| < 1e-6 m/s), `OutOfDomainReported`, `DeterministicBitwise`.

### Task 6.2: JSBSim model generation and runs
**Files:** Create `python/jsbsim_compare/{__init__.py, gen_xml.py, run_jsbsim.py, compare.py}`, `reference/jsbsim/{aircraft/f16_sl/f16_sl.xml, engine/dummy.xml, planet/flat_nonrot.xml}` (generated, committed), `scenarios/jsbsim/*.toml`, `python/tests/test_gen_xml.py`. Reuse knowledge from `spikes/jsbsim/` (FINDINGS.md, toy.xml) — re-implement in product code, don't import the spike.
- `gen_xml.py` reads `data/f16/f16_sl.toml` only. Body-axis X/Y/Z force and ROLL/PITCH/YAW moment functions reproducing §5.3 exactly (tables in deg breakpoints with `aero/alpha-deg`, `aero/beta-deg`; |β| fold via abs/sign functions; damping via `aero/bi2vel`, `aero/ci2vel`; aileron/rudder normalization; Cz's β² term with the 57.3 constant). CG transfer by geometry (AERORP at x_cg_ref·c̄, CG at x_cg·c̄, structural frame x aft, inches), NOT by the S&L formula. `<mass_balance negated_crossproduct_inertia="false">`; `emptywt` chosen so `inertia/mass-slugs` = 1/0.00157 (assert after load). Thrust: external body-x force at the CG = power-level interpolation of idle/mil/max tables at (altitude, Mach) with power P from a property set by the runner (constant at trim power; engine lag excluded per PRD §8.6). Dummy zero-power engine for FGTrim throttle; pitch trim → elevator pass-through with a known linear map; surfaces drivable directly. Rotor h_x = 0 on both sides for comparisons.
- Planet: spherical, `rotation_rate = 0`, radius R and GM such that g(0) = 9.80665 m/s²; `simulation/gravity-model = 0`; our side runs `gravity = {g0 = 9.80665, R = <same>}`, `atmosphere = "us76"`, `rotor = false`, `xcg = 0.35`.
- `run_jsbsim.py`: loads planet BEFORE model; sets integrators explicitly (rate: AB3 for rotational and translational; position: Trapezoidal / AB3 — record in metadata); dt = 1/2400 and 1/4800; writes CSV with the same column names as our `_fps/_deg` columns.
- [ ] `test_gen_xml.py`: every table value in the XML equals the TOML value (parse XML); mass/inertia properties after `load_model` equal the TOML values to 1e-9 rel; table-edge note: comparison scenarios assert in-domain.

### Task 6.3: Comparisons J1–J3 + discrepancy budget
**Files:** `python/jsbsim_compare/compare.py`, `python/tests/test_jsbsim_compare.py`, `reference/jsbsim/results/*.csv` (committed JSBSim outputs), report additions (`metrics.py`, `plots.py`, `docs/results.template.md`), CI job `jsbsim` in `.github/workflows/ci.yml`.
- J1: JSBSim initialized at our us76 trim state (V = 502 ft/s, h = 10,000 ft and h = 0 ft) → report `accelerations/udot-ft_sec2`, `wdot`, `qdot-rad_sec2` after `run_ic()`; criteria PRD §8.6 (initial targets; exceedance blocks until root-caused — do not loosen).
- J2: JSBSim's own trim vs ours: α, δe, throttle within 0.1 % or explained.
- J3: elevator doublet (±1°, 1–2 s / 2–3 s), aileron pulse (2°, 1–2 s), rudder doublet (±2°) from each simulator's own trim, 20 s, compare deviations from trim for ΔV, Δα, Δq, Δθ, Δh, Δβ, Δp, Δr, Δφ sampled at common times; metric = max |Δ_ours − Δ_jsb| / max|Δ_jsb| per signal; target ≤ 2 %.
- Discrepancy budget table: for each source (JSBSim dt halving; our dt halving; gravity model on/off in ours; atmosphere: ρ, a, T differences at trim altitude between our us76 and JSBSim's `atmosphere/*` properties; Earth curvature analytic estimate over the run; residual trim mismatch from J1; sampling alignment) quantify its effect on J3 metrics. Remaining unexplained difference stated plainly.
- CI job `jsbsim`: release build + `uv run python -m jsbsim_compare.compare --check` (regenerates XML, reruns both sides, asserts J1–J3 criteria and that regenerated JSBSim CSVs match committed ones within 1e-9 rel).
**Accept (M6):** PRD §13 M6.

## Milestone M7 — Actuators and autopilot

### Task 7.1: PID
**Files:** `src/flightsim/control/pid.hpp/.cpp`, `tests/unit/test_pid.cpp`.
`struct PidGains { double kp, ki, kd, tf /*derivative filter time constant, s*/, tt /*back-calculation tracking time, s*/, umin, umax; }; class Pid { public: explicit Pid(PidGains); double update(double setpoint, double measurement, double dt); void reset(double u0 = 0); bool antiwindup = true; };` Discrete form: P on error, D on measurement through first-order filter (Tustin or backward Euler — document), I with back-calculation `İ = ki·e + (u_sat − u)/tt`, output saturated. `reset(u0)` initializes the integrator so the first output equals u0 (bumpless engage from trim).
- [ ] Tests: `ProportionalOnly`, `IntegralRamp`, `DerivativeOnMeasurementNoKick` (setpoint step → no D spike), `SaturationClamps`, `BackCalculationLimitsWindup` (sustained saturation: integrator bounded vs unbounded with antiwindup=false), `BumplessReset`.

### Task 7.2: Autopilot + closed-loop runner
**Files:** `src/flightsim/control/autopilot.hpp/.cpp`; Modify `src/flightsim/sim/aircraft_runner.*`, `aircraft_scenario.*`; `scenarios/f16_ap_*.toml`; Test `tests/unit/test_autopilot.cpp`, `tests/validation/test_f16_autopilot.cpp`.
- Cascade per PRD §7.3: `h_cmd → PI → θ_cmd (±15°) → P → q_cmd → PI → δe_cmd`; `V_cmd → PI → throttle ∈ [0,1]`. Controller at fixed rate (default 50 Hz; `dt` must divide it — load error otherwise), sampled measurements (no noise), outputs held between samples; actuators enabled (lag/rate/position per data) and engine lag on. Engaged from trim with bumpless reset.
- Scenario `[autopilot] rate_hz, gains = "scenarios/autopilot_gains.toml", antiwindup = true`, `[[commands]] kind = "altitude"|"airspeed", t_s, value (+ units suffix)`, `[wind]`.
- Metrics computed in C++ (`control/metrics.hpp`: overshoot %, 2 % settling time, steady-state error over last 10 s, peak actuator rate/position) and written into metadata.json.
- [ ] Validation tests C1–C4 assert PRD §8.7 targets from the committed gains; C2 runs with and without anti-windup and asserts the anti-windup run has smaller overshoot.

### Task 7.3: Linear design + margins (Python)
**Files:** `python/control_design/{__init__.py, design.py}`, `python/tests/test_design.py`, `scenarios/autopilot_gains.toml` (generated, committed).
- From `flightsim linearize` JSON at the C-scenario trim: longitudinal model + first-order actuator lag + ZOH/Padé delay for the controller rate; design inner-to-outer loops (pitch-rate, pitch, altitude; airspeed) with python-control; report each loop's gain and phase margins broken one loop at a time (C5: GM ≥ 6 dB, PM ≥ 45°); write gains TOML; produce linear-predicted step responses for C1/C3 to overlay on the nonlinear runs.
- [ ] `test_design.py`: margins ≥ C5 for the committed gains; the committed gains file equals what `design.py` produces (reproducibility).

### Task 7.4: Report
Figures: C1–C3 responses with linear prediction overlay, C2 windup vs anti-windup, C4 gust/headwind rejection, actuator usage; results.md control table (targets vs measured).

**Accept (M7):** PRD §13 M7.

### (original M7 sketch, superseded) Interfaces: `control/pid.{hpp,cpp}`: `struct PidGains { double kp, ki, kd, tf, tt, umin, umax; }; class Pid { double update(double setpoint, double measurement, double dt); void reset(); bool antiwindup = true; }` (stateful by necessity — controllers are the one stateful component, outside derivative functions). `control/autopilot.{hpp,cpp}` cascade per §7.3; aircraft runner with controller scheduling (`sim/aircraft_runner.{hpp,cpp}`); `python/control_design/design.py` (margins C5). Scenarios C1–C4. **Accept:** PRD §13 M7.

## Milestone M8 — Stock JSBSim appendix + polish
J4 appendix, final results regeneration, README, decision log, limitations. **Accept:** PRD §13 M8.
