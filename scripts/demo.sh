#!/usr/bin/env bash
# One command: configure + build (release), run all tests, set up Python, run scenarios/studies/benchmarks,
# write out/index.html and docs/results.md. Exits non-zero on any failure.
set -euo pipefail
cd "$(dirname "$0")/.."
cmake --preset release
cmake --build --preset release -j
ctest --preset release --output-on-failure
cd python
uv sync --frozen
uv run python -m flightsim_analysis.run
uv run python -m flightsim_analysis.report
cd ..
echo "Report: $PWD/out/index.html"
