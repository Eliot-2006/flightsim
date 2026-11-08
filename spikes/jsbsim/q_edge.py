"""Table edge behaviour: 1-D table breakpoints alpha 0,0.1,0.2 -> 0,1,4."""
import common
from common import *
f = make(); init(f)
for a in (-0.1, 0.0, 0.15, 0.2, 0.3, 0.5):
    f["ic/alpha-rad"] = a; f.run_ic(); f.run()
    print(f"table edge: alpha={a:+.2f} -> {f['aero/edge-test']:.4f}")
