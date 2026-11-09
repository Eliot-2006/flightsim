# flightsim

A C++20 six-degree-of-freedom rigid-body simulator built to be checked, not just run. Every integrator is validated against closed-form solutions (point mass with drag, harmonic oscillator, torque-free rotation with Jacobi elliptic functions, the intermediate-axis flip), with measured convergence orders and conservation laws asserted in CI. An F-16 model, trim, linearization, JSBSim cross-validation and an altitude/airspeed autopilot are planned (see the roadmap below).

> This file is generated from `docs/README.template.md` by `python -m flightsim_analysis.report`; its numbers come from `out/metrics.json`.

## Headline results

{{orders.in_bounds}} of {{orders.asserted}} asserted convergence orders land inside their bounds: Euler, symplectic Euler and Lie–Euler measure {{orders.first.min:.3f}}–{{orders.first.max:.3f}} (bound [{{orders.first.lo}}, {{orders.first.hi}}]); RK4 and RKMK4 measure {{orders.fourth.min:.3f}}–{{orders.fourth.max:.3f}} (bound [{{orders.fourth.lo}}, {{orders.fourth.hi}}]).

![Convergence: error against step size for T2, T3, T4, R1, R2 and R3, with reference slopes](docs/img/convergence.png)

On the torque-free asymmetric body R3 at h = {{r3.h}} s over {{r3.t_final:.0f}} s, the maximum relative energy drift is {{r3.euler.dE:.1e}} for Euler, {{r3.rk4.dE:.1e}} for RK4 and {{r3.rkmk4.dE:.1e}} for RKMK4. Without renormalization the quaternion norm drifts {{r3.euler_norenorm.dq:.1e}} (Euler) and {{r3.rk4_norenorm.dq:.1e}} (RK4); RKMK4 needs no projection ({{r3.rkmk4.dq:.1e}}).

![Energy and angular-momentum drift on R3 per integrator](docs/img/drift.png)

Spin about the intermediate axis flips: the body rate q changes sign {{r4.sign_changes}} times in {{r4.t_final:.0f}} s while energy drifts {{r4.dE:.1e}} and angular momentum {{r4.dH:.1e}} (RK4, h = {{r4.h}} s).

![R4 Dzhanibekov flip: body rates and inertial angular momentum](docs/img/r4.png)

A rigid-body RK4 step with gravity and quadratic drag costs {{bench.step.rk4.ns:.0f}} ns ({{env.build_type}} build, {{env.compiler}}, {{bench.cpus}} logical CPUs (Google Benchmark reports no usable clock speed on this machine), single thread; Euler {{bench.step.euler.ns:.0f}} ns, RKMK4 {{bench.step.rkmk4.ns:.0f}} ns). The step loop performs zero heap allocations (asserted by `alloc_tests`).

## Quickstart

Requires CMake ≥ 3.28, a C++20 compiler, Ninja and [uv](https://docs.astral.sh/uv/). Dependencies (Eigen, GoogleTest, Google Benchmark, toml++) are fetched at configure time with pinned hashes.

```sh
scripts/demo.sh        # configure + build (release), ctest, uv sync --frozen, run, report
open out/index.html    # prints its path when done
```

Individual pieces:

```sh
cmake --preset release && cmake --build --preset release -j && ctest --preset release
build/release/flightsim run scenarios/rb_dzhanibekov.toml --out out/rb_dz        # trajectory.csv + metadata.json
build/release/flightsim converge scenarios/studies/core_convergence.toml --out out/convergence
build/release/bench_flightsim
cd python && uv sync --frozen && uv run pytest
```

`out/index.html` is a single static page with Plotly's JavaScript inlined (no network needed); it works in light and dark mode.

## Documents

- [`docs/PRD.md`](docs/PRD.md): the specification and validation strategy.
- [`docs/design.md`](docs/design.md): architecture, conventions with the test that pins each, derivations, decision log, limitations.
- [`docs/results.md`](docs/results.md): every measured number (generated).

## Roadmap

Done: rigid-body core, translational and rotational validation, five integrators (Euler, RK4, symplectic Euler, Lie–Euler, RKMK4), analysis tooling and report. Coming in M4–M7: F-16 model and consistency laws, trim and linearization, JSBSim cross-validation, actuators and autopilot.

MIT licensed.
