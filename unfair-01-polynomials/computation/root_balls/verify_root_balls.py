#!/usr/bin/env python3
"""Independent verifier for generate_root_balls.py output.

It recomputes every Rouché bound from exact H and rational JSON data.  A
stored success flag, stored slack, or generator record is never trusted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from fractions import Fraction
from pathlib import Path

from root_ball_common import (
    complex_from_json,
    complex_from_strings,
    disk_disjoint,
    disk_nonreal,
    fraction_string,
    lift_rouche_bounds,
    parse_fraction,
    polynomial_rouche_bounds,
)


H = [1, 1, 1, 1, 1, 1, 0, 0, 1]
L = 9


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def check_reported_bounds(group: dict[str, object], lhs: Fraction, rhs: Fraction) -> None:
    reported = group["rouche"]
    if parse_fraction(reported["lhs_upper"]) != lhs or parse_fraction(reported["rhs_lower"]) != rhs:
        raise AssertionError(f"stored Rouché bounds disagree for {group['id']}")
    if parse_fraction(reported["strict_slack_lower"]) != rhs - lhs:
        raise AssertionError(f"stored Rouché slack disagrees for {group['id']}")
    if not lhs < rhs:
        raise AssertionError(f"Rouché inequality fails for {group['id']}")


def get_disk(group: dict[str, object]):
    disk = group["representative"]
    radius = parse_fraction(disk["radius"])
    if radius <= 0:
        raise AssertionError(f"nonpositive disk radius: {group.get('id', '<unknown>')}")
    return complex_from_json(disk["center"]), radius


def all_disks(groups: list[dict[str, object]]):
    answer = []
    for group in groups:
        center, radius = get_disk(group)
        answer.append((group["id"] + ":representative", center, radius))
        answer.append((group["id"] + ":conjugate", (center[0], -center[1]), radius))
    return answer


def convolve(left: list[int], right: list[int]) -> list[int]:
    """Exact ascending-coefficient product, deliberately independent of input data."""
    answer = [0] * (len(left) + len(right) - 1)
    for i, coefficient in enumerate(left):
        for j, other in enumerate(right):
            answer[i + j] += coefficient * other
    return answer


def stretch(coefficients: list[int], exponent: int) -> list[int]:
    answer = [0] * (exponent * (len(coefficients) - 1) + 1)
    for index, coefficient in enumerate(coefficients):
        answer[exponent * index] = coefficient
    return answer


def coefficient_sha256(coefficients: list[int]) -> str:
    return hashlib.sha256(json.dumps(coefficients, separators=(",", ":")).encode()).hexdigest()


def require_fixed_r8_target(candidate: dict[str, object]) -> list[int]:
    """Bind the root-ball claim to the exact target R8(x)R8(x^9), not just a hash."""
    target = convolve(H, stretch(H, L))
    if candidate.get("H_coefficients_ascending") != H or candidate.get("L") != L:
        raise AssertionError("source candidate is not the fixed R8, L=9 instance")
    if candidate.get("degree") != len(target) - 1:
        raise AssertionError("source candidate has wrong target degree")
    if candidate.get("C_coefficients_ascending") != target:
        raise AssertionError("source candidate C is not the exact R8(x)R8(x^9) target")
    if candidate.get("C_sha256") != coefficient_sha256(target):
        raise AssertionError("source candidate C_sha256 does not match the exact fixed target")
    if any(coefficient not in (0, 1) for coefficient in target):
        raise AssertionError("fixed R8(x)R8(x^9) target unexpectedly is not binary")
    return target


def require_unique_source_coverage(groups: list[dict[str, object]], raw_groups: list[dict[str, object]]) -> dict[int, dict[str, object]]:
    """Ensure the 40 certificate groups cover each input quadratic group exactly once."""
    raw_by_index = {int(group["index"]): group for group in raw_groups}
    if len(raw_by_index) != len(raw_groups):
        raise AssertionError("source input repeats a quadratic group index")
    source_indices = [int(group["source_group_index"]) for group in groups]
    if len(set(source_indices)) != len(source_indices):
        raise AssertionError("duplicate source coverage in certificate groups")
    if set(source_indices) != set(raw_by_index):
        raise AssertionError("certificate source coverage is incomplete or has an unknown group")
    return raw_by_index


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    source = json.loads(args.input.read_text())
    candidate = source["candidates"][0]
    if certificate.get("format") != "R8-L9-rational-rouche-root-balls-v1":
        raise AssertionError("unexpected certificate format")
    if certificate.get("H_coefficients_ascending") != H or certificate.get("L") != L:
        raise AssertionError("wrong exact mother polynomial or L")
    if certificate.get("input_sha256") != sha(args.input) or certificate.get("input_candidate") != candidate["name"]:
        raise AssertionError("certificate does not bind to this source input")
    target = require_fixed_r8_target(candidate)
    scale = int(certificate["abs_bound_scale"])
    if scale <= 0:
        raise AssertionError("invalid modulus-bound scale")
    raw_groups = candidate["groups"]
    groups = certificate.get("groups")
    if not isinstance(groups, list) or len(groups) != 40:
        raise AssertionError("expected exactly 40 conjugate-pair groups")
    raw_by_index = require_unique_source_coverage(groups, raw_groups)
    ids = {group.get("id") for group in groups}
    expected_ids = {f"mother-{block}" for block in range(4)} | {f"lift-{block}-{lift}" for block in range(4) for lift in range(9)}
    if ids != expected_ids:
        raise AssertionError("group IDs are incomplete, duplicated, or unexpected")
    by_id = {group["id"]: group for group in groups}
    for block in range(4):
        group = by_id[f"mother-{block}"]
        raw = raw_by_index[int(group["source_group_index"])]
        if group["kind"] != "mother" or int(group["block"]) != block or raw["type"] != "original_quadratic":
            raise AssertionError(f"bad mother provenance for block {block}")
        center, radius = get_disk(group)
        if center != complex_from_strings(raw["parent_root_positive_imag"]):
            raise AssertionError(f"mother center was not copied exactly from source suggestion: block {block}")
        if group.get("conjugate_rule") != "(re,-im), same radius" or not disk_nonreal(center, radius):
            raise AssertionError(f"invalid mother conjugate/nonreal disk: block {block}")
        lhs, rhs = polynomial_rouche_bounds(H, center, radius, scale)
        check_reported_bounds(group, lhs, rhs)
    lift_parent_counts = {f"mother-{block}": 0 for block in range(4)}
    for block in range(4):
        parent = by_id[f"mother-{block}"]
        alpha_center, alpha_radius = get_disk(parent)
        for lift in range(9):
            group = by_id[f"lift-{block}-{lift}"]
            raw = raw_by_index[int(group["source_group_index"])]
            if (group["kind"] != "lift" or int(group["block"]) != block or int(group["lift_index"]) != lift or
                    group["parent_id"] != parent["id"] or raw["type"] != "lift_quadratic" or
                    int(raw["block"]) != block or int(raw["lift_index"]) != lift):
                raise AssertionError(f"bad lift provenance: block {block}, lift {lift}")
            center, radius = get_disk(group)
            if center != complex_from_strings(raw["lift_root_from_positive_imag_parent"]):
                raise AssertionError(f"lift center was not copied exactly from source suggestion: {group['id']}")
            if group.get("conjugate_rule") != "(re,-im), same radius" or not disk_nonreal(center, radius):
                raise AssertionError(f"invalid lift conjugate/nonreal disk: {group['id']}")
            if parse_fraction(group["rouche"]["alpha_radius_added"]) != alpha_radius:
                raise AssertionError(f"wrong mother-radius perturbation: {group['id']}")
            lhs, rhs = lift_rouche_bounds(alpha_center, alpha_radius, center, radius, scale)
            check_reported_bounds(group, lhs, rhs)
            lift_parent_counts[parent["id"]] += 1
    if set(lift_parent_counts.values()) != {9}:
        raise AssertionError("every certified mother root needs exactly nine lift disks")
    disks = all_disks(groups)
    if len(disks) != 80:
        raise AssertionError("expected 80 mother/lift disks including conjugates")
    closest_squared_slack = None
    for i, (name_a, center_a, radius_a) in enumerate(disks):
        if not disk_nonreal(center_a, radius_a):
            raise AssertionError(f"real-axis crossing: {name_a}")
        for name_b, center_b, radius_b in disks[:i]:
            if not disk_disjoint(center_a, radius_a, center_b, radius_b):
                raise AssertionError(f"overlap: {name_a}, {name_b}")
            squared_slack = ((center_a[0] - center_b[0]) ** 2 + (center_a[1] - center_b[1]) ** 2 -
                             (radius_a + radius_b) ** 2)
            if closest_squared_slack is None or squared_slack < closest_squared_slack:
                closest_squared_slack = squared_slack
    coverage = certificate.get("coverage", {})
    expected_coverage = {
        "mother_representatives": 4, "lift_representatives": 36, "conjugate_quadratic_groups": 40,
        "mother_root_disks": 8, "lift_root_disks": 72, "total_root_disks_for_H_times_H_x9": 80,
    }
    if any(coverage.get(key) != value for key, value in expected_coverage.items()) or coverage.get("pairwise_disjoint") is not True:
        raise AssertionError("incorrect coverage count")
    if parse_fraction(coverage["closest_squared_distance_minus_radius_sum_squared"]) != closest_squared_slack:
        raise AssertionError("stored closest disk-separation slack disagrees")
    # Mathematical coverage conclusion: Rouché gives one root in each of the
    # eight mother disks.  Degree(H)=8.  For the unique mother root alpha in
    # each disk, the nine lift inequalities give nine roots of z^9-alpha; the
    # degree is nine.  Pairwise disjointness yields the stated 8+72 coverage.
    record = {
        "verifier_sha256": sha(Path(__file__)), "common_sha256": sha(Path(__file__).with_name("root_ball_common.py")),
        "status": "verified", "certificate": str(args.certificate), "certificate_sha256": sha(args.certificate),
        "input": str(args.input), "input_sha256": sha(args.input), "groups": 40,
        "fixed_target": "R8(x)R8(x^9)", "fixed_target_degree": len(target) - 1,
        "fixed_target_coefficients_sha256": coefficient_sha256(target),
        "root_disks": 80, "all_rouche_inequalities_strict": True,
        "all_disks_nonreal_and_pairwise_disjoint": True,
        "coverage_conclusion": "8 simple roots of H and 72 simple roots of H(x^9), hence 80 separated roots for H(x)H(x^9)",
        "python": sys.version, "platform": platform.platform(), "float_arithmetic_used": False,
    }
    args.record.parent.mkdir(parents=True, exist_ok=True)
    atomic(args.record, record)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
