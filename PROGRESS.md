# PROGRESS

Living log for the 6-DOF flight dynamics simulator. Spec: `docs/PRD.md`. Plan: `docs/superpowers/plans/2026-11-08-flightsim.md`.

## Current milestone
M0 (foundations) — in progress. Parallel: M4.0 JSBSim spike, M4.1 F-16 data transcription.

## Orchestration
| Task | Track | Agent (model) | Branch / workspace | Status | Result |
|---|---|---|---|---|---|
| M0 foundations | A | Claude (Sonnet) | m0-foundations | running | — |
| 4.0 JSBSim spike | B | Claude (Sonnet) | spike-jsbsim | running | — |
| 4.1 F-16 data + reference values | C | Codex (gpt-6.1-sol) | f16-data | running | — |

Integration branch: `flight-dynamics-prd-81fca11f` → `main` via PR at verified milestone gates.

## Provider / quota notes
- Codex: GPT-6.1 Terra is not in the local model catalog; using `gpt-6.1-sol` (strongest available).
- Gemini CLI available (0.67.0) for optional second-opinion reviews. Antigravity not installed.
- No provider exposes remaining quota to this orchestrator; usage is kept conservative (≤ 3 concurrent workers, reviews only on math-heavy merges).

## Verification log
_(orchestrator-run checks only; worker claims are not recorded here until re-run)_

## Decisions
- 2026-11-08: Commits carry the owner's identity and no Co-Authored-By/“Generated with” trailers (owner's global CLAUDE.md overrides tool attribution defaults).
- 2026-11-08: Python pinned to 3.12 via uv (system Python 3.15 lacks jsbsim wheels risk).

## Problems and resolutions
_(none yet)_

## Remaining work
M0–M8 per plan.
