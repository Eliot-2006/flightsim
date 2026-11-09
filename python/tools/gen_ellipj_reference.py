"""Generate the Jacobi parameter-m reference with scipy.special.ellipj (A&S 16.4)."""

import csv
from pathlib import Path

import numpy as np
from scipy.special import ellipj


def main():
    path = Path(__file__).resolve().parents[2] / "tests/data/ellipj_reference.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["u", "m", "sn", "cn", "dn"])
        for m in (0.0, 0.1, 0.5, 0.9, 0.99, 0.999):
            for u in np.linspace(-5.0, 5.0, 101):
                sn, cn, dn, _ = ellipj(u, m)
                writer.writerow(format(float(v), ".17g") for v in (u, m, sn, cn, dn))
    print(f"wrote {path} (606 rows)")


if __name__ == "__main__":
    main()
