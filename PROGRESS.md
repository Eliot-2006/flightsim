# PROGRESS

Living log for the 6-DOF flight dynamics simulator. Spec: `docs/PRD.md`. Plan: `docs/superpowers/plans/2026-11-08-flightsim.md`.

## Current milestone
M8 (stock JSBSim appendix J4 + final polish) next. Done: M0–M7.

## Orchestration
| Task | Track | Agent (model) | Branch / workspace | Status | Result |
|---|---|---|---|---|---|
| M0 foundations | A | Claude (Sonnet) | m0-foundations | completed, merged | Presets dev/release/asan, pinned deps, CI matrix, uv env |
| M1 rigid-body core | A | Claude (Sonnet) | m1-core | completed, merged | 52 tests; T2 orders Euler 1.0015 / RK4 4.009; T5 catches transposed DCM and left-multiplication |
| M2 rotational + integrators | A | Codex (gpt-6.1-sol, high) | m2-rotational | completed, merged | 88 tests; 32/32 asserted orders in bounds (Euler ≈1.00–1.05, RK4/RKMK4 ≈4.00–4.01, symplectic 1.01) |
| M3 tooling, report, docs | A | Claude (Sonnet) | m3-tooling | completed, merged | demo.sh 28 s from clean clone; 113 tests; contrast ≥3:1 enforced by test |
| M4b wind, F-16 model, consistency | D | Codex (gpt-6.1-sol, high) | m4b-f16-model | completed, merged | A1 all 13 components ≤ 2e-4 rel (ṗ/ṙ gap = S&L rounded inertia coefficients); B1/B2 residuals 6e-11/3e-10 (order 4.0); S1 symmetrized 2e-16; aircraft RK4 order ≈2.2 |
| M6a JSBSim model generator | F | Claude (Sonnet) | m6a-jsbsim-model | completed, merged | Independent Python S&L evaluation vs JSBSim aero/thrust/mass ≤1e-9 rel; JSBSim trims at 0 and 10k ft; results at dt 1/2400 and 1/4800 committed |
| M6b runner + JSBSim J1–J3 + budget | F | Claude (Sonnet) | m6b-jsbsim-compare | completed, merged | J1 pass (ẇ offset = US76 ρ difference); J2 ≤ 1.7e-5 rel; J3 max 1.26 % (limit 2 %); remainder ≈ JSBSim's own dt-halving effect (~0.3 %) |
| M7a PID + linear design | G | Codex (gpt-6.1-sol, high) | m7a-pid-design | completed, merged | Linear C1: 4.3 % OS, 16.3 s settle, 3.0 ft/s; C3: 4.3 % OS, 4.2 ft; C5 GM ≥ 18 dB, PM ≥ 52° all loops (one airspeed→pitch cross-feed) |
| M7b autopilot closed loop | G | Claude (Sonnet, high) | m7b-autopilot | completed, merged | C1 4.47 % / 16.45 s / 3.19 ft/s (linear 4.30 % / 16.3 s); C2 AW 3.08 % / 21.2 s vs no-AW 4.61 % / 90.5 s; C3 4.22 %, 4.88 ft; C4 gust 17.1 ft; B1 closed-loop 1e-9 |
| M5 trim + linearization | E | Codex (gpt-6.1-sol, high) | m5-trim-linear | completed, merged | TR1: throttle +3.5e-5, δe +2e-5°, α +9e-7 rad (level_35); LN1 eigenvalue diffs ≤3.4e-4 1/s; LN2 ratio 96.3; LN3 ε=6e-6 flat |
| 4.2–4.3 interp, data loader, atmosphere | D | Claude (Sonnet) | m4a-data-atmos | completed, merged | 75 tests at merge; US76 vs NASA-TM-X-74335 ≤1e-4 |
| 4.0 JSBSim spike | B | Claude (Sonnet) | spike-jsbsim | completed, merged | Primary §8.6 approach feasible (jsbsim 1.3.3); see spikes/jsbsim/FINDINGS.md |
| 4.1 F-16 data + reference values | C | Codex (gpt-6.1-sol) | f16-data | completed, merged | 3 sources (AeroBench, Archimedes, UMN) + TP-1538 visual checks; 22 raw-A cells resolved; A1/TR1/LN1 single-lineage (Archimedes) |

Integration branch: `flight-dynamics-prd-81fca11f` → `main`. The owner opened PR #2 (branch → main) and merged PR #1 themselves, so the orchestrator keeps the branch green and leaves PR merges to the owner.

## Provider / quota notes
- Codex: GPT-6.1 Terra is not in the local model catalog; using `gpt-6.1-sol` (strongest available).
- Gemini CLI available (0.67.0) for optional second-opinion reviews. Antigravity not installed.
- Codex shows its 5h limit in the TUI (30% left at 16:04, reset 16:59; 50% left at 17:49, resets 22:00); Claude/Gemini expose none.
- No other provider exposes remaining quota to this orchestrator; usage is kept conservative (≤ 3 concurrent workers, reviews only on math-heavy merges).

## Verification log
_(orchestrator-run checks only; worker claims are not recorded here until re-run)_
- 2026-11-08 M0: clean `release` configure+build+ctest in m0 worktree: 3/3 pass; `uv run pytest`: 1 passed.
- 2026-11-08 M1: release build, ctest 52/52 (41 unit, 11 validation); CLI run writes CSV+JSON; reviewed quat/EOM/RK4 by hand. CI (Linux gcc/clang dev+release, asan, format, python, bench) green on M0 commit.
- 2026-11-08 M2+M4a merged: integration release build ctest 111/111 (75 unit, 36 validation); `flightsim converge` orders.json re-generated and inspected.
- 2026-11-08 M7 merged: release ctest 202/202 (AppleClang + GCC 13); pytest 40 passed; demo.sh exit 0; jsbsim compare --check passed; C2 figure inspected.
- 2026-11-08 CI adacc6d: all jobs green except sanitizers pending at check time (cross-platform regression tolerances fixed).
- 2026-11-08 M7a merged: release ctest 172/172 (AppleClang + GCC 13); pytest 40 passed.
- 2026-11-08 CI fixes 43e551d: JSBSim CSV regression check now uses per-unit physical floors + ψ wrap (Linux libm round-off in unexcited lateral channels); alloc test skipped under sanitizers (operator new replacement vs ASan); clang-format 18 pinned via uvx.
- 2026-11-08 M6 merged: release ctest 165/165 (AppleClang + GCC 13 container); `uv run pytest` 36 passed; `jsbsim_compare.compare --check` CHECK PASSED.
- 2026-11-08 M5+M6a merged: release ctest 159/159 (AppleClang), 159/159 (GCC 13 container); `scripts/demo.sh` exit 0; `uv run pytest` 23 passed.
- 2026-11-08 M4: integration release ctest 139/139 (AppleClang) and 139/139 (GCC 13 container).
- 2026-11-08 M3: fresh clone `scripts/demo.sh` exit 0 in 28.4 s, 113/113 tests; figures inspected (convergence slopes, R4 flip with constant H_N).
- 2026-11-08 Local GCC 13 container check (`gcc_check.sh`) added to the merge routine; integration branch passes.
- 2026-11-08 CI: gcc-13 failed on -Wdangling-else in a test (fixed d8849c9); macos-14 runner label unsupported → macos-15.
- 2026-11-08 4.1: re-ran f16_data_diff.py → PASS, 6 documented resolutions, 0 unresolved.
- 2026-11-08 4.0: re-ran `spikes/jsbsim/q2_3_4.py`; integrator/planet outputs match FINDINGS.md.

## Decisions
- 2026-11-08: Commits carry the owner's identity and no Co-Authored-By/“Generated with” trailers (owner's global CLAUDE.md overrides tool attribution defaults).
- 2026-11-08: Python pinned to 3.12 via uv.
- 2026-11-08: Autopilot operating point 502 ft/s at 10,000 ft (SlAdc, xcg 0.35, actuators + engine lag on, 50 Hz) — keeps gust/step dips away from the ground stop.
- 2026-11-08: J3 runs at 10,000 ft (from sea level the maneuvers leave the thrust table's lower edge where S&L extrapolates and JSBSim clamps).
- 2026-11-08: LN1 will use the xcg = 0.30 published eigenvalues (only published linear data found); TR1 uses xcg 0.30/0.35/0.38 level trims; all labelled single-lineage per PRD §8.5.
- 2026-11-08: PRD §10.3 amended: Euler on R1 has attitude order 2 (exact rotation axis), asserted on norm error instead.
- 2026-11-08: PRD §8.6 amended per spike (no-RK4 JSBSim reference via two small steps; J4 trim-only).

## Problems and resolutions
- AppleClang 15 ASan aborts at startup on this Mac → local asan runs use Homebrew LLVM (`CC/CXX=$(brew --prefix llvm)/bin/clang[++]`); CI uses clang-18 on Linux.
- `gh` CLI binary is x86 and won't run on this arm64 Mac → CI status read via public GitHub API; PRs created via web/API later.
- System /usr/local pythons are x86 → uv-managed Python 3.12.
- Codex workers can stall on a first-run 'Trust this folder?' prompt → orchestrator checks each Codex terminal right after launch and answers it.
- GCC-only -Wdangling-else failures from unbraced gtest macros (twice) → braced; worker rules updated; local GCC 13 container check before every merge.
- CMake 3.25.1 mis-parses `SYSTEM` in FetchContent_Declare → minimum raised to 3.28 (tested in CI and container).
- M4b BLOCKED (correctly): planned 60 s open-loop fixture left the α table domain at 25 s; same-method 1e-3 reference rule infeasible for the non-smooth aircraft → PRD amended (20 s in-domain run; 10× reference rule for aircraft; 60 s closed-loop budget check in M7).
- M5 BLOCKED on LN2 breakpoint monitor (β crosses 0 via rotor coupling) → criterion clarified (inactive/smooth crossings excluded with in-test proofs); ratio bound unchanged.
- M7 review: first C2 anti-windup run never settled (tt = 1/ωc dragged the integrator) → tt = Ti per PI loop (Åström–Murray rule) via design.py; C2 criterion tightened in PRD.
- M7a BLOCKED twice: gains file picked up as a report scenario (moved under scenarios/control); first design missed linear C1/C3 → redesigned in design.py (not by trial on the nonlinear model).
- Codex 5h limit hit during M4b → kept model, resumed after reset.
- JSBSim has no RK4; stock f16 FCS overrides direct surface writes → PRD §8.6 amendment (J4 trim-only).

## Remaining work
M0–M8 per plan.
