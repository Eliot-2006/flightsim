"""Run the flightsim CLI over every scenario and study, plus the benchmarks. Python never integrates anything."""

import argparse
import os
import json
import re
import tomllib
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
    sh(bench, f"--benchmark_min_time={args.bench_min_time}", f"--benchmark_out={out / 'bench.json'}",
       "--benchmark_out_format=json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
