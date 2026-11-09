#!/usr/bin/env python3
"""Compare F-16 numeric transcriptions; run --self-test for regression checks."""
import argparse
import math
from pathlib import Path
import sys
import tomllib
import unittest


def differences(left, right, path=""):
    if isinstance(left, dict) and isinstance(right, dict):
        result = []
        for key in sorted((left.keys() | right.keys()) - {"source"}):
            name = f"{path}.{key}" if path else key
            if key not in left or key not in right:
                result.append((name, left.get(key, "<missing>"), right.get(key, "<missing>")))
            else:
                result.extend(differences(left[key], right[key], name))
        return result
    if isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            return [(path + ".length", len(left), len(right))]
        return [item for i, (a, b) in enumerate(zip(left, right))
                for item in differences(a, b, f"{path}[{i}]")]
    if type(left) in (int, float) and type(right) in (int, float):
        # Task 4.1 specifies relative 1e-9, with no absolute tolerance.
        equal = (math.isfinite(left) and math.isfinite(right)
                 and math.isclose(left, right, rel_tol=1e-9, abs_tol=0.0))
    else:
        equal = type(left) is type(right) and left == right
    return [] if equal else [(path, left, right)]


def resolve(mismatches, resolutions):
    if not isinstance(resolutions, list):
        raise ValueError("resolutions must be an array of tables")
    pending = {path: (a, b) for path, a, b in mismatches}
    for record in resolutions:
        if (not isinstance(record, dict)
                or any(not isinstance(record.get(k), str) or not record[k].strip()
                       for k in ("path", "reason", "source"))
                or any(type(record.get(k)) not in (int, float)
                       or not math.isfinite(record[k]) for k in ("left", "right", "chosen"))):
            raise ValueError("resolution requires path, reason, source and finite numeric values")
        path = record["path"]
        if path not in pending:
            raise ValueError(f"stale, duplicate, or undocumented resolution: {path}")
        a, b = pending[path]
        if (type(a) not in (int, float) or type(b) not in (int, float)
                or record["left"] != a or record["right"] != b
                or record["chosen"] != a):
            raise ValueError(f"resolution does not pin both values and select canonical left: {path}")
        del pending[path]
    return [(path, a, b) for path, (a, b) in pending.items()]


class ComparisonTests(unittest.TestCase):
    def test_reports_changed_cell_and_missing_data(self):
        self.assertEqual(differences({"cx": [[1.0, 2.0]]}, {"cx": [[1.0, 3.0]]}),
                         [("cx[0][1]", 2.0, 3.0)])
        self.assertTrue(differences({"cx": [1]}, {}))
        self.assertTrue(differences({"cx": [1]}, {"cx": [1, 2]}))

    def test_relative_tolerance_has_no_absolute_floor(self):
        # Task 4.1 requires relative 1e-9; zero must agree exactly.
        self.assertFalse(differences(1.0, 1.0 + 5e-10))
        self.assertTrue(differences(1.0, 1.0 + 2e-9))
        self.assertTrue(differences(0.0, 1e-15))

    def test_only_citations_are_ignored(self):
        self.assertFalse(differences({"source": "A", "units": "deg"},
                                     {"source": "B", "units": "deg"}))
        self.assertTrue(differences({"units": "deg"}, {"units": "rad"}))
        self.assertTrue(differences(True, 1))
        self.assertTrue(differences(float("nan"), float("nan")))
        self.assertTrue(differences(float("inf"), float("inf")))

    def test_resolution_pins_both_values_and_requires_evidence(self):
        mismatch = [("cz.values[3]", -0.415, -0.416)]
        record = {"path": "cz.values[3]", "left": -0.415, "right": -0.416,
                  "chosen": -0.415, "reason": "Verified reproduction", "source": "citation"}
        self.assertEqual(resolve(mismatch, [record]), [])
        for field, value in [("left", -0.414), ("right", -0.417), ("reason", ""),
                             ("source", ""), ("chosen", 0.0)]:
            with self.subTest(field=field), self.assertRaises(ValueError):
                resolve(mismatch, [{**record, field: value}])
        with self.assertRaises(ValueError):
            resolve([], [record])
        with self.assertRaises(ValueError):
            resolve(mismatch, [record, record])


    def test_malformed_or_nonfinite_resolutions_cannot_hide_errors(self):
        record = {"path": "x", "left": 1.0, "right": 2.0, "chosen": 1.0,
                  "reason": "Evidence", "source": "citation"}
        for records in [None, {}, [None], [{**record, "left": True}],
                        [{**record, "reason": 3}]]:
            with self.subTest(records=records), self.assertRaises(ValueError):
                resolve([("x", 1.0, 2.0)], records)
        with self.assertRaises(ValueError):
            resolve([("x", float("inf"), 2.0)],
                    [{**record, "left": float("inf"), "chosen": float("inf")}])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", type=Path, nargs="?")
    parser.add_argument("right", type=Path, nargs="?")
    parser.add_argument("--resolutions", type=Path,
                        help="defaults to resolutions.toml beside LEFT, if present")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return not unittest.TextTestRunner().run(
            unittest.defaultTestLoader.loadTestsFromTestCase(ComparisonTests)).wasSuccessful()
    if args.left is None or args.right is None:
        parser.error("LEFT and RIGHT are required unless --self-test is used")
    try:
        with args.left.open("rb") as stream:
            left = tomllib.load(stream)
        with args.right.open("rb") as stream:
            right = tomllib.load(stream)
        records = []
        resolution_path = args.resolutions or args.left.with_name("resolutions.toml")
        if args.resolutions or resolution_path.exists():
            with resolution_path.open("rb") as stream:
                records = tomllib.load(stream)["resolution"]
        mismatches = differences(left, right)
        unresolved = resolve(mismatches, records)
        for path, a, b in mismatches:
            status = "MISMATCH" if (path, a, b) in unresolved else "RESOLVED"
            print(f"{status} {path}: {a!r} != {b!r}")
        print(f"{'FAIL' if unresolved else 'PASS'}: {len(mismatches)} differences, "
              f"{len(records)} documented resolutions, {len(unresolved)} unresolved")
        return bool(unresolved)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
