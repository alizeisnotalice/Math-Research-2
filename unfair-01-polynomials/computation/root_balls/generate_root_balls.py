#!/usr/bin/env python3
"""Generate an exact-rational Rouché disk certificate for R8--L9 roots.

The decimal roots in the input are used solely as rational centre suggestions.
All acceptance decisions below are made with Fraction arithmetic and integer
square roots; no float or numerical library is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

from root_ball_common import (
    ComplexQ,
    abs_bounds,
    complex_from_strings,
    complex_to_json,
    disk_disjoint,
    disk_nonreal,
    fraction_string,
    lift_rouche_bounds,
    parse_fraction,
    polynomial_rouche_bounds,
)
from fractions import Fraction


H = [1, 1, 1, 1, 1, 1, 0, 0, 1]  # ascending exact coefficients
ABS_SCALE = 10**160
MOTHER_RADIUS = Fraction(1, 10**70)
LIFT_RADIUS = Fraction(1, 10**65)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def disk_json(center: ComplexQ, radius: Fraction) -> dict[str, object]:
    return {"center": complex_to_json(center), "radius": fraction_string(radius)}


def add_mother(raw: dict, block: int) -> dict[str, object]:
    center = complex_from_strings(raw["parent_root_positive_imag"])
    lhs, rhs = polynomial_rouche_bounds(H, center, MOTHER_RADIUS, ABS_SCALE)
    if not lhs < rhs:
        raise RuntimeError(f"mother block {block} fails Rouché: {lhs} >= {rhs}")
    if not disk_nonreal(center, MOTHER_RADIUS):
        raise RuntimeError(f"mother block {block} disk meets real axis")
    return {
        "id": f"mother-{block}", "kind": "mother", "block": block,
        "source_group_index": int(raw["index"]), "representative": disk_json(center, MOTHER_RADIUS),
        "conjugate_rule": "(re,-im), same radius", "rouche": {
            "polynomial": "H", "lhs_upper": fraction_string(lhs),
            "rhs_lower": fraction_string(rhs), "strict_slack_lower": fraction_string(rhs - lhs),
        },
    }


def add_lift(raw: dict, block: int, mother: dict[str, object]) -> dict[str, object]:
    center = complex_from_strings(raw["lift_root_from_positive_imag_parent"])
    mother_center = mother["representative"]["center"]
    alpha_center = (parse_fraction(mother_center["re"]), parse_fraction(mother_center["im"]))
    alpha_radius = parse_fraction(mother["representative"]["radius"])
    lhs, rhs = lift_rouche_bounds(alpha_center, alpha_radius, center, LIFT_RADIUS, ABS_SCALE)
    if not lhs < rhs:
        raise RuntimeError(f"lift block {block}, index {raw['lift_index']} fails Rouché: {lhs} >= {rhs}")
    if not disk_nonreal(center, LIFT_RADIUS):
        raise RuntimeError(f"lift block {block}, index {raw['lift_index']} disk meets real axis")
    return {
        "id": f"lift-{block}-{int(raw['lift_index'])}", "kind": "lift", "block": block,
        "lift_index": int(raw["lift_index"]), "source_group_index": int(raw["index"]),
        "parent_id": mother["id"], "representative": disk_json(center, LIFT_RADIUS),
        "conjugate_rule": "(re,-im), same radius", "rouche": {
            "polynomial": "z^9-alpha for alpha in parent mother disk",
            "alpha_radius_added": fraction_string(alpha_radius),
            "lhs_upper": fraction_string(lhs), "rhs_lower": fraction_string(rhs),
            "strict_slack_lower": fraction_string(rhs - lhs),
        },
    }


def all_disks(groups: list[dict[str, object]]) -> list[tuple[str, ComplexQ, Fraction]]:
    answer = []
    for group in groups:
        disk = group["representative"]
        center_data = disk["center"]
        center = parse_fraction(center_data["re"]), parse_fraction(center_data["im"])
        radius = parse_fraction(disk["radius"])
        answer.append((group["id"] + ":representative", center, radius))
        answer.append((group["id"] + ":conjugate", (center[0], -center[1]), radius))
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    args = parser.parse_args()
    if not Path("/Volumes/data").is_mount():
        raise SystemExit("/Volumes/data is not mounted; refusing to write elsewhere")
    raw_input = json.loads(args.input.read_text())
    candidate = raw_input["candidates"][0]
    if candidate["H_coefficients_ascending"] != H or candidate["L"] != 9:
        raise SystemExit("expected exact R8 mother H and L=9 input")
    raw_groups = candidate["groups"]
    mothers_raw = {int(group["block"]): group for group in raw_groups if group["type"] == "original_quadratic"}
    if set(mothers_raw) != {0, 1, 2, 3}:
        raise SystemExit("expected four mother groups indexed by blocks 0..3")
    groups = [add_mother(mothers_raw[block], block) for block in range(4)]
    mother_by_id = {group["id"]: group for group in groups}
    lifts = [group for group in raw_groups if group["type"] == "lift_quadratic"]
    if len(lifts) != 36:
        raise SystemExit("expected 36 lift groups")
    for raw in sorted(lifts, key=lambda item: (int(item["block"]), int(item["lift_index"]))):
        block = int(raw["block"])
        groups.append(add_lift(raw, block, mother_by_id[f"mother-{block}"]))
    disks = all_disks(groups)
    if len(disks) != 80:
        raise AssertionError("expected 80 root disks")
    closest_squared_slack = None
    for i, (name_a, center_a, radius_a) in enumerate(disks):
        if not disk_nonreal(center_a, radius_a):
            raise AssertionError(f"{name_a} meets real axis")
        for name_b, center_b, radius_b in disks[:i]:
            if not disk_disjoint(center_a, radius_a, center_b, radius_b):
                raise AssertionError(f"overlapping disks: {name_a}, {name_b}")
            squared_slack = ((center_a[0] - center_b[0]) ** 2 + (center_a[1] - center_b[1]) ** 2 -
                             (radius_a + radius_b) ** 2)
            if closest_squared_slack is None or squared_slack < closest_squared_slack:
                closest_squared_slack = squared_slack
    certificate = {
        "format": "R8-L9-rational-rouche-root-balls-v1",
        "claim_boundary": "Exact rational Rouché disk certificate only; it does not certify coefficient screening or branch-and-bound exclusion.",
        "H_coefficients_ascending": H, "L": 9, "input_sha256": sha(args.input),
        "input_candidate": candidate["name"], "abs_bound_scale": str(ABS_SCALE),
        "groups": groups,
        "coverage": {
            "mother_representatives": 4, "lift_representatives": 36,
            "conjugate_quadratic_groups": 40, "mother_root_disks": 8,
            "lift_root_disks": 72, "total_root_disks_for_H_times_H_x9": 80,
            "pairwise_disjoint": True, "closest_squared_distance_minus_radius_sum_squared": fraction_string(closest_squared_slack),
            "argument": "Four disjoint mother disks plus conjugates give eight H roots. For each, nine disjoint lift disks use the z^9-alpha Rouché inequality; degree counts cover H(x^9). Cross-disjointness gives 80 simple roots for H(x)H(x^9).",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.record.parent.mkdir(parents=True, exist_ok=True)
    atomic(args.output, certificate)
    record = {
        "status": "generated", "certificate": str(args.output), "certificate_sha256": sha(args.output),
        "input": str(args.input), "input_sha256": sha(args.input), "groups": len(groups), "disks": len(disks),
        "python": sys.version, "platform": platform.platform(), "float_arithmetic_used": False,
    }
    atomic(args.record, record)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
