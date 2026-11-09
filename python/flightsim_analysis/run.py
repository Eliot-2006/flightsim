"""Run the flightsim CLI over every scenario and study, plus the benchmarks. Python never integrates anything."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bin", default=os.environ.get("FLIGHTSIM_BIN", ROOT / "build" / "release" / "flightsim"))
    ap.add_argument("--out", default=ROOT / "out")
    ap.add_argument("--bench-min-time", default="0.2s")
    args = ap.parse_args(argv)
    exe, out = Path(args.bin), Path(args.out)

    def sh(*cmd):
        print("+", *cmd, flush=True)
        subprocess.run([str(c) for c in cmd], check=True)

    for scenario in sorted((ROOT / "scenarios").glob("*.toml")):
        sh(exe, "run", scenario, "--out", out / "runs" / scenario.stem)
    for study in sorted((ROOT / "scenarios" / "studies").glob("*.toml")):
        sh(exe, "converge", study, "--out", out / "studies" / study.stem)
    bench = exe.parent / "bench_flightsim"
    sh(bench, f"--benchmark_min_time={args.bench_min_time}", f"--benchmark_out={out / 'bench.json'}",
       "--benchmark_out_format=json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
