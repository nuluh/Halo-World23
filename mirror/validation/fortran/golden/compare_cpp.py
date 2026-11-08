#!/usr/bin/env python3
"""Compare modern C++ stage probe values with frozen Fortran oracle fixtures."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
GOLDEN = Path(__file__).resolve().parent
CATALOG = GOLDEN / "catalog.json"

STAGE_ORDER = [
    ("basis", lambda k: k.startswith("basis.")),
    ("mass_inertia", lambda k: k.startswith("tmass.") or k.startswith("rmas_node")),
    ("element_axial", lambda k: k.startswith("axial_extension.")),
    ("element_torsion", lambda k: k.startswith("torsion.")),
    ("element_bend_y", lambda k: k.startswith("bend_y_node1.")),
    ("element_bend_z", lambda k: k.startswith("bend_z_node1.")),
    ("element_combined", lambda k: k.startswith("combined.")),
    ("element_rigid_body", lambda k: k.startswith("rigid_z.")),
    (
        "single_integration_step",
        lambda k: k.startswith(("dprev_isolated.", "dnext_isolated.",
                                "velocity_isolated.", "accel_isolated.",
                                "increment_isolated.")),
    ),
    ("continuous_beam", lambda k: k.startswith("continuous_beam.")),
    ("planar_portal", lambda k: k.startswith("planar_portal.")),
    ("spatial_frame", lambda k: k.startswith("spatial_frame.")),
    ("short_transient", lambda k: k.startswith("short_transient.")),
]


class CompareError(RuntimeError):
    pass


def load_csv(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if not row:
                continue
            if len(row) != 2:
                raise CompareError(f"{path}: expected key,value CSV")
            key = row[0].strip()
            if key in values:
                raise CompareError(f"{path}: duplicate key {key}")
            values[key] = float(row[1].strip())
    return values


def unit_for(key: str, fixture: dict) -> str:
    matches = [
        rule["unit_class"]
        for rule in fixture["unit_rules"]
        if re.fullmatch(rule["pattern"], key)
    ]
    if len(matches) != 1:
        raise CompareError(
            f"{fixture['id']}: {key!r} matched {len(matches)} unit rules: {matches}"
        )
    return matches[0]


def stage_for(key: str) -> str:
    for name, predicate in STAGE_ORDER:
        if predicate(key):
            return name
    raise CompareError(f"probe emitted unclassified stage key: {key}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("probe_csv", type=Path)
    args = parser.parse_args()

    with CATALOG.open(encoding="utf-8") as handle:
        catalog = json.load(handle)

    expected: dict[str, tuple[float, dict]] = {}
    for fixture in catalog["fixtures"]:
        if fixture["id"] not in {
            "3dframe-basic-basis-mass-load",
            "3dframe-element-force",
            "3dframe-central-step",
            "3dframe-progressive-frames-transient",
        }:
            continue
        values = load_csv(ROOT / fixture["output"]["path"])
        for key, value in values.items():
            expected[key] = (value, fixture)

    actual = load_csv(args.probe_csv)
    policy = catalog["tolerance_policy"]
    rtol = policy["relative_tolerance"]

    stage_stats: dict[str, dict[str, float | int]] = {}
    first_failure: tuple[str, str, float, float, float] | None = None

    for key, got in actual.items():
        if key not in expected:
            raise CompareError(f"probe key has no frozen oracle value: {key}")
        want, fixture = expected[key]
        unit = unit_for(key, fixture)
        definition = policy["unit_classes"][unit]
        if definition["comparison"] != "numeric":
            continue
        atol = definition["absolute_floor"]
        limit = atol + rtol * abs(want)
        err = abs(got - want)
        stage = stage_for(key)
        stats = stage_stats.setdefault(stage, {"count": 0, "max_abs_error": 0.0})
        stats["count"] = int(stats["count"]) + 1
        stats["max_abs_error"] = max(float(stats["max_abs_error"]), err)
        if err > limit and first_failure is None:
            first_failure = (stage, key, got, want, limit)

    for name, _ in STAGE_ORDER:
        if name in stage_stats:
            stats = stage_stats[name]
            print(
                f"[{name}] {stats['count']} values; "
                f"max_abs_error={float(stats['max_abs_error']):.17e}"
            )

    if first_failure is not None:
        stage, key, got, want, limit = first_failure
        print(
            "C++/Fortran parity FAILED at earliest stage "
            f"{stage}: {key}: got={got:.17e}, expected={want:.17e}, "
            f"allowed_abs_error={limit:.17e}",
            file=sys.stderr,
        )
        return 1

    print(
        f"C++/Fortran parity OK: {len(actual)} modern stage values "
        "match frozen oracle evidence"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError, CompareError) as exc:
        print(f"C++/Fortran parity ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
