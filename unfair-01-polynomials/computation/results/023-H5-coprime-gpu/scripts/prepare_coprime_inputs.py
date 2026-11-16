#!/usr/bin/env python3
"""Prepare numerical root blocks for H(x^a)H(x^b), gcd(a,b)=1.

The output is a new schema: its root blocks contain only a- and b-lift
quadratics, and every option records selected global group indices explicitly.
It is numerical preparation, not an algebraic certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from math import comb, gcd
from pathlib import Path

import mpmath as mp
import numpy as np

H = [1, 1, 1, 0, 0, 1, 1]
HERE = Path(__file__).resolve().parent
mp.mp.dps = 100


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integer_conv(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] += x * y
    return answer


def sparse_substitute(coefficients: list[int], exponent: int) -> list[int]:
    answer = [0] * ((len(coefficients) - 1) * exponent + 1)
    answer[::exponent] = coefficients
    return answer


def fmt(value: mp.mpf) -> str:
    return mp.nstr(value, 80)


def pair(poly: list[mp.mpf]) -> dict:
    values = np.asarray([float(value) for value in poly], dtype=np.float64)
    high = values.astype(np.float32)
    return {"hi": high.tolist(), "lo": (values - high.astype(np.float64)).astype(np.float32).tolist(),
            "float64": values.tolist()}


def conv_mp(left: list[mp.mpf], right: list[mp.mpf]) -> list[mp.mpf]:
    answer = [mp.mpf("0")] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] += x * y
    return answer


def quadratic(root: mp.mpc) -> list[mp.mpf]:
    return [abs(root) ** 2, -2 * mp.re(root), mp.mpf(1)]


def positive_imaginary_mother_roots() -> tuple[list[mp.mpc], str]:
    roots = mp.polyroots(list(reversed(H)), maxsteps=1000, error=False)
    positive = sorted((root for root in roots if mp.im(root) > 0), key=lambda z: (float(mp.re(z)), float(mp.im(z))))
    assert len(positive) == 3
    residual = max(abs(sum(coefficient * root ** i for i, coefficient in enumerate(H))) /
                   sum(coefficient * abs(root) ** i for i, coefficient in enumerate(H)) for root in roots)
    return positive, fmt(residual)


def local_options(a: int, b: int, u: int, v: int) -> list[tuple[int, tuple[int, ...], tuple[int, ...]]]:
    """Selections ordered by their combined local bit mask."""
    records = []
    for choose_a in itertools.combinations(range(a), u):
        for choose_b in itertools.combinations(range(b), v):
            positions = tuple(choose_a) + tuple(a + j for j in choose_b)
            mask = sum(1 << position for position in positions)
            records.append((mask, tuple(choose_a), tuple(choose_b)))
    return sorted(records, key=lambda item: item[0])


def product(indices: list[int], groups: list[dict]) -> list[mp.mpf]:
    answer = [mp.mpf(1)]
    for index in indices:
        answer = conv_mp(answer, groups[index]["coefficients_mp"],)
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=int, required=True)
    parser.add_argument("--b", type=int, default=7)
    parser.add_argument("--out-input", type=Path, required=True)
    parser.add_argument("--out-options", type=Path, required=True)
    args = parser.parse_args()
    a, b = args.a, args.b
    if not (1 < a < b and gcd(a, b) == 1):
        raise SystemExit("require 1<a<b and gcd(a,b)=1")
    if not Path("/Volumes/data").is_mount():
        raise SystemExit("/Volumes/data is not mounted; refusing to write elsewhere")
    roots, mother_residual = positive_imaginary_mother_roots()
    source_groups: list[dict] = []
    blocks: list[dict] = []
    for block_number, alpha in enumerate(roots):
        a_indices, b_indices = [], []
        for family, count, saved in (("a_lift", a, a_indices), ("b_lift", b, b_indices)):
            radius = abs(alpha) ** (mp.mpf(1) / count)
            theta = mp.arg(alpha)
            for lift_index in range(count):
                root = radius * mp.e ** (1j * (theta + 2 * mp.pi * lift_index) / count)
                index = len(source_groups)
                source_groups.append({"index": index, "block": block_number, "family": family,
                                      "lift_index": lift_index, "coefficients_mp": quadratic(root),
                                      "coefficients_ascending": [float(value) for value in quadratic(root)],
                                      "root_positive_imag_parent": [fmt(mp.re(alpha)), fmt(mp.im(alpha))]})
                saved.append(index)
        blocks.append({"block": block_number,
                       "parent_root_positive_imag": [fmt(mp.re(alpha)), fmt(mp.im(alpha))],
                       "a_lift_groups": a_indices, "b_lift_groups": b_indices,
                       "all_group_indices": a_indices + b_indices})
    C = integer_conv(sparse_substitute(H, a), sparse_substitute(H, b))
    if set(C) - {0, 1}:
        raise SystemExit("H(x^a)H(x^b) is not binary")
    terminal = a * b
    by_E: dict[int, tuple[int, int]] = {}
    for u in range(a + 1):
        for v in range(b + 1):
            E = b * u + a * v
            if 0 < E < terminal:
                if E in by_E:
                    raise SystemExit(f"nonunique internal balance E={E}")
                by_E[E] = (u, v)
    expected_E = sorted(by_E)
    cache_blocks = []
    for block in blocks:
        group_map = {group["index"]: group for group in source_groups}
        options: dict[str, list[dict]] = {}
        for E in expected_E:
            u, v = by_E[E]
            rows = []
            for local_mask, selected_a, selected_b in local_options(a, b, u, v):
                selected = [block["a_lift_groups"][j] for j in selected_a] + [block["b_lift_groups"][j] for j in selected_b]
                complement = [index for index in block["all_group_indices"] if index not in set(selected)]
                rows.append({"mode": "mixed", "u": u, "v": v, "selected_local_mask_decimal": str(local_mask),
                             "selected_group_indices": selected,
                             "selected_local_pairs": ([{"family": "a_lift", "lift_index": j, "group_index": block["a_lift_groups"][j]} for j in selected_a] +
                                                      [{"family": "b_lift", "lift_index": j, "group_index": block["b_lift_groups"][j]} for j in selected_b]),
                             "Ablock": pair(product(selected, group_map)), "Bblock": pair(product(complement, group_map))})
            assert len(rows) == comb(a, u) * comb(b, v)
            options[str(E)] = rows
        all_a = list(block["a_lift_groups"])
        all_b = list(block["b_lift_groups"])
        options[str(terminal)] = [
            {"mode": "all_a", "u": a, "v": 0, "selected_group_indices": all_a,
             "selected_local_pairs": [{"family": "a_lift", "lift_index": j, "group_index": all_a[j]} for j in range(a)],
             "Ablock": pair(product(all_a, group_map)), "Bblock": pair(product(all_b, group_map))},
            {"mode": "all_b", "u": 0, "v": b, "selected_group_indices": all_b,
             "selected_local_pairs": [{"family": "b_lift", "lift_index": j, "group_index": all_b[j]} for j in range(b)],
             "Ablock": pair(product(all_b, group_map)), "Bblock": pair(product(all_a, group_map))},
        ]
        cache_blocks.append({"block": block["block"], "a_lift_groups": block["a_lift_groups"],
                             "b_lift_groups": block["b_lift_groups"], "all_group_indices": block["all_group_indices"],
                             "options": options})
    # Strip mpmath-only fields from the portable raw input after all options exist.
    groups = [{key: value for key, value in group.items() if key != "coefficients_mp"} for group in source_groups]
    franel = lambda n: sum(comb(n, k) ** 3 for k in range(n + 1))
    total = franel(a) * franel(b) // 2 + 2
    raw = {"name": f"coprime_H6_a{a}_b{b}", "kind": "coprime_substitution",
           "H_coefficients_ascending": H, "a": a, "b": b, "degree": len(C) - 1,
           "C_coefficients_ascending": C, "C_sha256": hashlib.sha256(json.dumps(C, separators=(",", ":")).encode()).hexdigest(),
           "root_precision_decimal_digits": mp.mp.dps, "mother_root_max_relative_residual": mother_residual,
           "parent_block_count": len(blocks), "total_groups": len(groups), "groups": groups, "blocks": blocks,
           "balance_terminal_E": terminal,
           "balance_equation": "E=b*u+a*v; all parent blocks have identical E",
           "E_options": [{"E": E, "u": by_E[E][0], "v": by_E[E][1], "per_block_options": comb(a, by_E[E][0]) * comb(b, by_E[E][1])} for E in expected_E] +
                        [{"E": terminal, "special_modes": ["all_a", "all_b"], "per_block_options": 2, "first_block_fixed_mode": "all_a"}],
           "balanced_enumeration": {"constant_and_full_excluded": True, "complement_canonical": True,
                                     "Franel_a": franel(a), "Franel_b": franel(b), "nonconstant_complement_representatives": total},
           "status": "numerical_prepared_not_algebraically_certified"}
    args.out_input.parent.mkdir(parents=True, exist_ok=True)
    args.out_input.write_text(json.dumps({"candidates": [raw]}, indent=2) + "\n")
    cache = {"format": "coprime-balanced-ds-block-options-v1", "candidate": raw["name"], "kind": raw["kind"],
             "a": a, "b": b, "balance_terminal_E": terminal, "mp_dps": mp.mp.dps, "block_count": len(cache_blocks),
             "blocks": cache_blocks, "input_sha256": sha(args.out_input), "source_sha256": sha(Path(__file__)),
             "mixed_radix_order": "block order; first block radix is least significant; terminal fixes its all_a option",
             "claim_boundary": "Numerical root blocks only; root balancing is a search restriction, not an exact sign certificate."}
    args.out_options.write_text(json.dumps(cache, indent=2) + "\n")
    print(json.dumps({"candidate": raw["name"], "a": a, "b": b, "E_options": [row["E"] for row in raw["E_options"]],
                      "total": total, "input_sha256": sha(args.out_input), "options_sha256": sha(args.out_options)}, indent=2))


if __name__ == "__main__":
    main()
