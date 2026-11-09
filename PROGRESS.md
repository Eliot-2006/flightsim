# PROGRESS

Living log for the 6-DOF flight dynamics simulator. Spec: `docs/PRD.md`. Plan: `docs/superpowers/plans/2026-11-08-flightsim.md`.

## Current milestone
M2 complete (MVP). Next: M3 (tooling/docs/bench) and M4b (F-16 model, Tasks 4.4–4.6) in parallel. Done: M0, M1, M2, M4.0, M4.1, M4.2–4.3.

## Orchestration
| Task | Track | Agent (model) | Branch / workspace | Status | Result |
|---|---|---|---|---|---|
| M0 foundations | A | Claude (Sonnet) | m0-foundations | completed, merged | Presets dev/release/asan, pinned deps, CI matrix, uv env |
| M1 rigid-body core | A | Claude (Sonnet) | m1-core | completed, merged | 52 tests; T2 orders Euler 1.0015 / RK4 4.009; T5 catches transposed DCM and left-multiplication |
| M2 rotational + integrators | A | Codex (gpt-6.1-sol, high) | m2-rotational | completed, merged | 88 tests; 32/32 asserted orders in bounds (Euler ≈1.00–1.05, RK4/RKMK4 ≈4.00–4.01, symplectic 1.01) |
| 4.2–4.3 interp, data loader, atmosphere | D | Claude (Sonnet) | m4a-data-atmos | completed, merged | 75 tests at merge; US76 vs NASA-TM-X-74335 ≤1e-4 |
| 4.0 JSBSim spike | B | Claude (Sonnet) | spike-jsbsim | completed, merged | Primary §8.6 approach feasible (jsbsim 1.3.3); see spikes/jsbsim/FINDINGS.md |
| 4.1 F-16 data + reference values | C | Codex (gpt-6.1-sol) | f16-data | completed, merged | 3 sources (AeroBench, Archimedes, UMN) + TP-1538 visual checks; 22 raw-A cells resolved; A1/TR1/LN1 single-lineage (Archimedes) |

Integration branch: `flight-dynamics-prd-81fca11f` → `main` via PR at verified milestone gates.

## Provider / quota notes
- Codex: GPT-6.1 Terra is not in the local model catalog; using `gpt-6.1-sol` (strongest available).
- Gemini CLI available (0.67.0) for optional second-opinion reviews. Antigravity not installed.
- Codex shows its 5h limit in the TUI (30% left at 16:04, resets 16:59); Claude/Gemini expose none.
- No other provider exposes remaining quota to this orchestrator; usage is kept conservative (≤ 3 concurrent workers, reviews only on math-heavy merges).

## Verification log
_(orchestrator-run checks only; worker claims are not recorded here until re-run)_
- 2026-11-08 M0: clean `release` configure+build+ctest in m0 worktree: 3/3 pass; `uv run pytest`: 1 passed.
- 2026-11-08 M1: release build, ctest 52/52 (41 unit, 11 validation); CLI run writes CSV+JSON; reviewed quat/EOM/RK4 by hand. CI (Linux gcc/clang dev+release, asan, format, python, bench) green on M0 commit.
- 2026-11-08 M2+M4a merged: integration release build ctest 111/111 (75 unit, 36 validation); `flightsim converge` orders.json re-generated and inspected.
- 2026-11-08 CI: gcc-13 failed on -Wdangling-else in a test (fixed d8849c9); macos-14 runner label unsupported → macos-15.
- 2026-11-08 4.1: re-ran f16_data_diff.py → PASS, 6 documented resolutions, 0 unresolved.
- 2026-11-08 4.0: re-ran `spikes/jsbsim/q2_3_4.py`; integrator/planet outputs match FINDINGS.md.

## Decisions
- 2026-11-08: Commits carry the owner's identity and no Co-Authored-By/“Generated with” trailers (owner's global CLAUDE.md overrides tool attribution defaults).
- 2026-11-08: Python pinned to 3.12 via uv.
- 2026-11-08: LN1 will use the xcg = 0.30 published eigenvalues (only published linear data found); TR1 uses xcg 0.30/0.35/0.38 level trims; all labelled single-lineage per PRD §8.5.
- 2026-11-08: PRD §10.3 amended: Euler on R1 has attitude order 2 (exact rotation axis), asserted on norm error instead.
- 2026-11-08: PRD §8.6 amended per spike (no-RK4 JSBSim reference via two small steps; J4 trim-only).

## Problems and resolutions
- AppleClang 15 ASan aborts at startup on this Mac → local asan runs use Homebrew LLVM (`CC/CXX=$(brew --prefix llvm)/bin/clang[++]`); CI uses clang-18 on Linux.
- `gh` CLI binary is x86 and won't run on this arm64 Mac → CI status read via public GitHub API; PRs created via web/API later.
- System /usr/local pythons are x86 → uv-managed Python 3.12.
- Codex workers can stall on a first-run 'Trust this folder?' prompt → orchestrator checks each Codex terminal right after launch and answers it.
- JSBSim has no RK4; stock f16 FCS overrides direct surface writes → PRD §8.6 amendment (J4 trim-only).

## Remaining work
M0–M8 per plan.
