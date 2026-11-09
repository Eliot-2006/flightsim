# Project record

Final record of how flightsim was built and verified. Spec: [`docs/PRD.md`](docs/PRD.md) (with every amendment). Plan: [`docs/superpowers/plans/2026-11-08-flightsim.md`](docs/superpowers/plans/2026-11-08-flightsim.md). Results: [`docs/results.md`](docs/results.md). Live report: <https://eliot-2006.github.io/flightsim/>.

## Status

Complete. All PRD milestones (M0–M8) meet their acceptance criteria. An independent review found no physics or maths defects, and its five code findings are fixed. Every CI job is green: GCC 13, Clang 18 and AppleClang builds in debug and release, ASan/UBSan, format, Python, JSBSim comparison, benchmarks, coverage (gated at 85 %), and a fresh-clone demo on Linux. The report is published to GitHub Pages from `main`.

Not built (stretch goals in the PRD, documented as future work): heading hold, an LQR comparison, Dryden turbulence, C¹ table interpolation, and an in-browser (WASM) simulator.

## Milestones

| Milestone | Outcome |
|---|---|
| M0 Foundations | CMake presets (dev/release/asan), pinned dependencies by SHA-256, CI matrix, uv-managed Python 3.12 |
| M1 Rigid-body core | Quaternion/frame conventions, each with a test; 13-state EOM; Euler and RK4; TOML scenarios; CLI `run`; translational analytic cases T0–T5 |
| M2 Rotational validation | Symplectic Euler, Lie–Euler, RKMK4; torque-free references R1–R5 (Jacobi elliptic); 32/32 asserted convergence orders inside their bounds; spin-stability linearization |
| M3 Tooling | Benchmarks, zero-allocation test, Plotly report with generated results.md and README, one-command `scripts/demo.sh` |
| M4 F-16 model | S&L data transcribed and cross-checked against three sources; interpolation, atmospheres, engine, actuators, wind; check case A1 ≤ 2e-4; energy and angular-impulse balances 6e-11 / 3e-10 |
| M5 Trim and linearization | Trim within 0.028 % of published values; eigenvalues within 0.28 %; LN2 quadratic scaling ratio 96 |
| M6 JSBSim | Generated JSBSim aircraft from the same data on a matched planet; J1–J3 pass (worst 1.26 % of peak against 2 %); discrepancy budget |
| M7 Autopilot | Linear design with margins ≥ 18 dB / ≥ 53°; nonlinear C1–C4 pass; anti-windup cuts C2 settling from 90 s to 21 s |
| M8 Polish | Stock-JSBSim F-16 appendix (J4, trim only), F-16 benchmarks, coverage and demo CI, final documentation |
| Release pass | Recruiter-facing README and report (results at a glance, chapters, 3D flight path from simulator output), GitHub Pages |

## How it was built

The work was orchestrated across agents in isolated Git worktrees, and every branch was re-verified by the orchestrator before merging: a release and debug build, the full test suite, a GCC 13 container build, and pytest. Claude (Sonnet) built M0, M1, M3, the JSBSim tooling, the closed-loop autopilot and the polish. Codex (gpt-6.1-sol) built the numerically heavy parts: M2, the F-16 model, trim and linearization, and the control design. Claude (Opus) did the independent review; the Gemini CLI could not be used (the client is no longer supported).

## Decisions

Each one is recorded with its rationale in the decision log in [`docs/design.md`](docs/design.md).

- Flat, non-rotating Earth; SI internally; velocity in NED; scalar-first Hamilton quaternion body→NED.
- No book access, so the F-16 tables come from open reproductions cross-checked between sources, with NASA TP-1538 as the primary citation. Published reference values are labelled single-lineage throughout.
- JSBSim comparison against a JSBSim aircraft generated from the same data file, on a non-rotating planet with matched gravity; J3 is flown at 10,000 ft so both simulators stay inside the thrust tables.
- Autopilot operating point: 502 ft/s at 10,000 ft. Gains are designed only in `python/control_design/design.py`, never tuned on the nonlinear model.
- PRD amendments, each with its evidence: Euler's order-2 attitude convergence on R1; the reference rule for the non-smooth aircraft model; the 20 s open-loop B1/B2 runs plus a 60 s closed-loop budget check; the LN2 breakpoint criterion; the tightened C2 criterion; and JSBSim's lack of RK4.

## Problems found and fixed

- A 60 s open-loop maneuver left the α table domain, so B1/B2 run for 20 s open-loop and 60 s closed-loop.
- Back-calculation anti-windup at tt = 1/ωc never settled C2. Fixed with tt = Ti, derived in the design script.
- The P-only PID back-calculation left a bias; anti-windup tracked the commanded elevator instead of the rate-limited one. Both fixed from the independent review, with regression tests.
- Cross-platform floating-point differences: regression checks against committed JSBSim CSVs and gains compare numerically with physical floors, never bitwise.
- Toolchain: GCC-only `-Wdangling-else` in tests; CMake 3.25 mis-parsing `SYSTEM` (minimum raised to 3.28); ASan versus a replaced global `operator new`; a GitHub release-asset 504 (googletest now fetched from the tag archive, same SHA-256).
