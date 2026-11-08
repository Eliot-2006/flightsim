# PROGRESS

Living log for the 6-DOF flight dynamics simulator. Spec: `docs/PRD.md`. Plan: `docs/superpowers/plans/2026-11-08-flightsim.md`.

## Current milestone
M1 (rigid-body core) — in progress. Parallel: M4.1 F-16 data transcription. Done: M0, M4.0 spike.

## Orchestration
| Task | Track | Agent (model) | Branch / workspace | Status | Result |
|---|---|---|---|---|---|
| M0 foundations | A | Claude (Sonnet) | m0-foundations | completed, merged | Presets dev/release/asan, pinned deps, CI matrix, uv env |
| M1 rigid-body core | A | Claude (Sonnet) | m1-core | running | — |
| 4.0 JSBSim spike | B | Claude (Sonnet) | spike-jsbsim | completed, merged | Primary §8.6 approach feasible (jsbsim 1.3.3); see spikes/jsbsim/FINDINGS.md |
| 4.1 F-16 data + reference values | C | Codex (gpt-6.1-sol) | f16-data | running | — |

Integration branch: `flight-dynamics-prd-81fca11f` → `main` via PR at verified milestone gates.

## Provider / quota notes
- Codex: GPT-6.1 Terra is not in the local model catalog; using `gpt-6.1-sol` (strongest available).
- Gemini CLI available (0.67.0) for optional second-opinion reviews. Antigravity not installed.
- No provider exposes remaining quota to this orchestrator; usage is kept conservative (≤ 3 concurrent workers, reviews only on math-heavy merges).

## Verification log
_(orchestrator-run checks only; worker claims are not recorded here until re-run)_
- 2026-11-08 M0: clean `release` configure+build+ctest in m0 worktree: 3/3 pass; `uv run pytest`: 1 passed.
- 2026-11-08 4.0: re-ran `spikes/jsbsim/q2_3_4.py`; integrator/planet outputs match FINDINGS.md.

## Decisions
- 2026-11-08: Commits carry the owner's identity and no Co-Authored-By/“Generated with” trailers (owner's global CLAUDE.md overrides tool attribution defaults).
- 2026-11-08: Python pinned to 3.12 via uv.
- 2026-11-08: PRD §8.6 amended per spike (no-RK4 JSBSim reference via two small steps; J4 trim-only).

## Problems and resolutions
- AppleClang 15 ASan aborts at startup on this Mac → local asan runs use Homebrew LLVM (`CC/CXX=$(brew --prefix llvm)/bin/clang[++]`); CI uses clang-18 on Linux.
- `gh` CLI binary is x86 and won't run on this arm64 Mac → CI status read via public GitHub API; PRs created via web/API later.
- System /usr/local pythons are x86 → uv-managed Python 3.12.
- JSBSim has no RK4; stock f16 FCS overrides direct surface writes → PRD §8.6 amendment (J4 trim-only).

## Remaining work
M0–M8 per plan.
