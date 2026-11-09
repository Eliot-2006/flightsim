"""Run the flightsim CLI over every scenario and study, plus the benchmarks. Python never integrates anything."""

import argparse
import os
import json
import platform
import re
import tomllib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# Benchmark settings (reported in docs/results.md): the median of BENCH_REPETITIONS repetitions, each at least
# --bench-min-time long; per-repetition iteration counts are chosen by Google Benchmark.
BENCH_REPETITIONS = 7
BUILD_FLAGS = ("-O", "-DNDEBUG", "-std=", "-ffp-contract", "-march", "-arch")  # compile-command tokens worth reporting


# (test binary, output name, gtest filter) of the tests whose printed numbers metrics.checks() reads
CHECKS = (("validation_tests", "validation", "F16Consistency.*:F16AutopilotB1.*:F16Trim.TR2Stationarity"),
          ("unit_tests", "unit", "F16Model.A1CheckCase"))


def cpu_model() -> str:
    """CPU brand string: sysctl on macOS, /proc/cpuinfo on Linux (Google Benchmark's MHz is unusable on Apple silicon)."""
    if platform.system() == "Darwin":
        return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True,
                              check=True).stdout.strip()
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name"):
            return line.split(":", 1)[1].strip()
    return platform.processor() or "unknown"


def environment(build_dir: Path, bench_min_time: str) -> dict:
    """Hardware, OS, compiler, build type and flags of the build that produced the benchmark numbers."""
    cache = (build_dir / "CMakeCache.txt").read_text()
    var = lambda name: re.search(rf"^{name}:[A-Z]+=(.*)$", cache, re.MULTILINE).group(1)
    compiler = var("CMAKE_CXX_COMPILER")
    version = subprocess.run([compiler, "--version"], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    command = next(e["command"] for e in json.loads((build_dir / "compile_commands.json").read_text())
                   if e["file"].endswith("bench/bench_f16.cpp"))
    tokens = command.split()
    flags = [" ".join(tokens[i:i + 2]) if t == "-arch" else t for i, t in enumerate(tokens) if t.startswith(BUILD_FLAGS)]
    return {"cpu": cpu_model(), "logical_cpus": os.cpu_count(), "os": platform.platform(), "compiler": version,
            "build_type": var("CMAKE_BUILD_TYPE"), "flags": " ".join(flags),
            "bench_min_time": bench_min_time, "bench_repetitions": BENCH_REPETITIONS}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bin", default=os.environ.get("FLIGHTSIM_BIN", ROOT / "build" / "release" / "flightsim"))
    ap.add_argument("--out", default=ROOT / "out")
    ap.add_argument("--bench-min-time", default="0.1s")
    args = ap.parse_args(argv)
    exe, out = Path(args.bin), Path(args.out)

    def sh(*cmd):
        print("+", *cmd, flush=True)
        subprocess.run([str(c) for c in cmd], check=True)

    for scenario in sorted((ROOT / "scenarios").glob("*.toml")):
        config = tomllib.loads(scenario.read_text())
        if "autopilot" in config:  # closed-loop C1-C4: run (trim/linearize reject actuator-enabled models)
            sh(exe, "run", scenario, "--out", out / "control" / scenario.stem)
        elif "aircraft" in config:
            tag = f"level_{round(config['aircraft']['xcg'] * 100):02d}"
            sh(exe, "trim", scenario, "--out", out / "aircraft" / f"trim_{tag}.json")
            sh(exe, "linearize", scenario, "--out", out / "aircraft" / f"linear_{tag}.json")
        else:
            sh(exe, "run", scenario, "--out", out / "runs" / scenario.stem)
    # The third published CG is the same trim specification with xcg=.38; no extra physics in Python.
    generated = out / "aircraft" / "f16_trim_502_xcg038.toml"
    config_text = (ROOT / "scenarios" / "f16_trim_502_xcg035.toml").read_text()
    config_text = re.sub(r"^xcg\s*=.*$", "xcg = 0.38", config_text, flags=re.MULTILINE)
    config_text = re.sub(r"^data\s*=.*$", "data = " + json.dumps(str(ROOT / "data/f16/f16_sl.toml")),
                         config_text, flags=re.MULTILINE)
    generated.write_text(config_text)
    sh(exe, "trim", generated, "--out", out / "aircraft" / "trim_level_38.json")
    # C5 margins and the linear step predictions come from the design script (drives the same CLI's linearizer).
    from control_design import design

    design.generate(out / "control_design")
    for study in sorted((ROOT / "scenarios" / "studies").glob("*.toml")):
        sh(exe, "converge", study, "--out", out / "studies" / study.stem)
    # JSBSim matched comparison (J1-J3 + budget): our runs via the CLI, JSBSim side from reference/jsbsim.
    from jsbsim_compare import compare

    compare.main(["--out", str(out / "jsbsim"), "--bin", str(exe)])
    bench = exe.parent / "bench_flightsim"
    # Numbers the validation tests print (A1, B1, B2, S1, W1, TR2, aircraft convergence, closed-loop B1): re-run just
    # those tests and keep their stdout for metrics.py. The assertions themselves are enforced by ctest.
    (out / "checks").mkdir(parents=True, exist_ok=True)
    for binary, name, gtest_filter in CHECKS:
        with open(out / "checks" / f"{name}.txt", "w") as f:
            print("+", binary, flush=True)
            subprocess.run([str(exe.parent / binary), f"--gtest_filter={gtest_filter}"], check=True, stdout=f)
    # Stock JSBSim f16 trim (J4 appendix): trim only, a different model from the matched one above.
    from jsbsim_compare import stock_f16

    stock_f16.main(["--out", str(out / "jsbsim"), "--bin", str(exe)])
    sh(bench, f"--benchmark_min_time={args.bench_min_time}", f"--benchmark_repetitions={BENCH_REPETITIONS}",
       "--benchmark_report_aggregates_only=true", f"--benchmark_out={out / 'bench.json'}", "--benchmark_out_format=json")
    (out / "env.json").write_text(json.dumps(environment(exe.parent, args.bench_min_time), indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
