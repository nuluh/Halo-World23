#!/usr/bin/env python3
"""Validate golden Fortran-oracle metadata, hashes, units, and tolerance evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
GOLDEN = Path(__file__).resolve().parent
CATALOG_PATH = GOLDEN / "catalog.json"
SCHEMA_PATH = GOLDEN / "schema.json"


class VerificationError(RuntimeError):
    pass


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def check_hash(path_string: str, expected: str, label: str) -> None:
    path = ROOT / path_string
    require(path.is_file(), f"{label}: missing file: {path_string}")
    actual = sha256(path)
    require(actual == expected, f"{label}: SHA-256 mismatch for {path_string}: {actual} != {expected}")


def read_fixture_csv(path_string: str) -> list[tuple[str, float]]:
    path = ROOT / path_string
    rows: list[tuple[str, float]] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle), start=1):
            if not row:
                continue
            require(len(row) == 2, f"{path_string}:{line_number}: expected exactly two CSV fields")
            key = row[0].strip()
            require(key and key not in seen, f"{path_string}:{line_number}: duplicate/empty key {key!r}")
            seen.add(key)
            try:
                value = float(row[1].strip())
            except ValueError as exc:
                raise VerificationError(f"{path_string}:{line_number}: non-numeric value") from exc
            require(math.isfinite(value), f"{path_string}:{line_number}: non-finite value")
            rows.append((key, value))
    require(rows, f"{path_string}: empty fixture")
    return rows


def classify(key: str, rules: list[dict], fixture_id: str) -> str:
    matches = [rule["unit_class"] for rule in rules if re.fullmatch(rule["pattern"], key)]
    require(
        len(matches) == 1,
        f"{fixture_id}: key {key!r} must match exactly one unit rule, matched {matches}",
    )
    return matches[0]


def verify_catalog(catalog: dict) -> None:
    require(catalog.get("schema_version") == 1, "unsupported catalog schema_version")
    env = catalog["reference_environment"]
    env_id = env["id"]

    policy = catalog["tolerance_policy"]
    multiplier = policy["relative_multiplier_epsilon"]
    expected_rtol = multiplier * sys.float_info.epsilon
    require(
        math.isclose(policy["relative_tolerance"], expected_rtol, rel_tol=0.0, abs_tol=1e-30),
        "relative_tolerance does not equal multiplier * DBL_EPSILON",
    )
    unit_classes = policy["unit_classes"]

    observed: dict[str, float] = {name: 0.0 for name in unit_classes}
    fixture_ids: set[str] = set()

    for fixture in catalog["fixtures"]:
        fid = fixture["id"]
        require(fid not in fixture_ids, f"duplicate fixture id: {fid}")
        fixture_ids.add(fid)
        require(
            fixture["reference_environment_id"] == env_id,
            f"{fid}: unknown reference environment",
        )

        for source in fixture["legacy_sources"]:
            check_hash(source["path"], source["sha256"], fid)

        generation = fixture["generation"]
        require(
            isinstance(generation.get("compiler"), str) and bool(generation["compiler"]),
            f"{fid}: missing generation compiler",
        )
        require(
            isinstance(generation.get("compiler_version"), str)
            and bool(generation["compiler_version"]),
            f"{fid}: missing generation compiler version",
        )
        check_hash(generation["driver"], generation["driver_sha256"], fid)
        check_hash(generation["runner"], generation["runner_sha256"], fid)
        check_hash(fixture["output"]["path"], fixture["output"]["sha256"], fid)

        rows = read_fixture_csv(fixture["output"]["path"])
        for key, value in rows:
            unit_class = classify(key, fixture["unit_rules"], fid)
            require(unit_class in unit_classes, f"{fid}: undefined unit class {unit_class!r}")
            observed[unit_class] = max(observed[unit_class], abs(value))

    for unit_class, definition in unit_classes.items():
        max_expected = definition["max_abs_expected"]
        require(
            math.isclose(observed[unit_class], max_expected, rel_tol=1e-15, abs_tol=1e-30),
            f"{unit_class}: max_abs_expected stale: catalog={max_expected}, observed={observed[unit_class]}",
        )
        expected_floor = policy["relative_tolerance"] * max(1.0, max_expected)
        require(
            math.isclose(definition["absolute_floor"], expected_floor, rel_tol=1e-15, abs_tol=1e-30),
            f"{unit_class}: absolute_floor does not follow policy formula",
        )

    evidence = {entry["category"]: entry for entry in catalog["evidence_groups"]}
    required_categories = {
        "geometry_local_basis",
        "lumped_mass_inertia",
        "element_forces",
        "global_force_assembly",
        "single_integration_step",
        "short_transient_history",
        "longer_regression_runs",
    }
    require(set(evidence) == required_categories, "tolerance evidence categories are incomplete")

    for category, entry in evidence.items():
        if entry["status"] == "measured":
            require(
                entry["observed_max_abs_error"] is not None,
                f"{category}: measured evidence requires an observed error",
            )
            require(
                entry["planned_issue"] is None,
                f"{category}: measured evidence must not retain a planned issue",
            )
        elif entry["status"] == "pending_measurement":
            require(
                entry["observed_max_abs_error"] is None,
                f"{category}: pending evidence must not invent an observed error",
            )
            require(
                isinstance(entry["planned_issue"], int),
                f"{category}: pending evidence requires a downstream issue",
            )
        else:
            raise VerificationError(f"{category}: unsupported evidence status {entry['status']!r}")

    # Presence/readability check for the published JSON Schema. Full JSON-Schema
    # validation is intentionally not required so this verifier stays stdlib-only.
    schema = load_json(SCHEMA_PATH)
    require(schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema", "unexpected schema dialect")


def regenerate(catalog: dict) -> None:
    for fixture in catalog["fixtures"]:
        command = fixture["generation"]["command"]
        runner = fixture["generation"]["runner"]
        expected = f"bash {runner}"
        require(command == expected, f"{fixture['id']}: regeneration command must be exactly {expected!r}")
        print(f"[regenerate] {fixture['id']}: {command}", flush=True)
        subprocess.run(["bash", str(ROOT / runner)], cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--regenerate",
        action="store_true",
        help="compile/run preserved Fortran drivers and byte-compare regenerated snapshots",
    )
    args = parser.parse_args()

    try:
        catalog = load_json(CATALOG_PATH)
        verify_catalog(catalog)
        if args.regenerate:
            regenerate(catalog)
            verify_catalog(catalog)
    except (VerificationError, OSError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(f"golden-fixture verification FAILED: {exc}", file=sys.stderr)
        return 1

    print(
        f"golden-fixture verification OK: {len(catalog['fixtures'])} fixtures, "
        f"{len(catalog['tolerance_policy']['unit_classes'])} unit classes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
