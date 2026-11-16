#!/usr/bin/env python3
"""Independent exact/layout verifier for the coprime substitution cache."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from math import comb
from pathlib import Path

import numpy as np


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def conv(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] += x * y
    return answer


def substitute(h: list[int], q: int) -> list[int]:
    result = [0] * ((len(h) - 1) * q + 1)
    result[::q] = h
    return result


def unpair(row: dict) -> np.ndarray:
    return np.asarray(row["hi"], dtype=np.float64) + np.asarray(row["lo"], dtype=np.float64)


def franel(n: int) -> int:
    return sum(comb(n, k) ** 3 for k in range(n + 1))


def local_records(a: int, b: int, u: int, v: int) -> list[tuple[int, tuple[int, ...], tuple[int, ...]]]:
    answer = []
    for ca in itertools.combinations(range(a), u):
        for cb in itertools.combinations(range(b), v):
            positions = tuple(ca) + tuple(a + index for index in cb)
            answer.append((sum(1 << index for index in positions), tuple(ca), tuple(cb)))
    return sorted(answer, key=lambda row: row[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--options", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads(args.input.read_text())["candidates"][0]
    cache = json.loads(args.options.read_text())
    h = raw["H_coefficients_ascending"]
    a, b = int(raw["a"]), int(raw["b"])
    terminal = a * b
    assert raw["kind"] == cache["kind"] == "coprime_substitution"
    assert cache["candidate"] == raw["name"] and (cache["a"], cache["b"], cache["balance_terminal_E"]) == (a, b, terminal)
    C = conv(substitute(h, a), substitute(h, b))
    assert C == raw["C_coefficients_ascending"] and set(C) <= {0, 1}
    assert len(raw["blocks"]) == len(cache["blocks"]) == raw["parent_block_count"] == 3
    assert len(raw["groups"]) == raw["total_groups"] == 3 * (a + b)
    group_by_index = {int(group["index"]): group for group in raw["groups"]}
    balance: dict[int, tuple[int, int]] = {}
    for u in range(a + 1):
        for v in range(b + 1):
            E = b * u + a * v
            if 0 < E < terminal:
                assert E not in balance
                balance[E] = (u, v)
    expected_E = sorted(balance) + [terminal]
    assert [item["E"] for item in raw["E_options"]] == expected_E
    layouts: dict[int, list[int]] = {}
    for raw_block, cached_block in zip(raw["blocks"], cache["blocks"]):
        assert raw_block["block"] == cached_block["block"]
        ai = list(raw_block["a_lift_groups"])
        bi = list(raw_block["b_lift_groups"])
        assert cached_block["a_lift_groups"] == ai and cached_block["b_lift_groups"] == bi
        assert raw_block["all_group_indices"] == ai + bi == cached_block["all_group_indices"]
        for j, index in enumerate(ai):
            assert group_by_index[index]["block"] == raw_block["block"] and group_by_index[index]["family"] == "a_lift" and group_by_index[index]["lift_index"] == j
        for j, index in enumerate(bi):
            assert group_by_index[index]["block"] == raw_block["block"] and group_by_index[index]["family"] == "b_lift" and group_by_index[index]["lift_index"] == j
        assert sorted(map(int, cached_block["options"])) == expected_E
        for E in sorted(balance):
            u, v = balance[E]
            rows = cached_block["options"][str(E)]
            expected = local_records(a, b, u, v)
            assert len(rows) == len(expected) == comb(a, u) * comb(b, v)
            layouts.setdefault(E, []).append(len(rows))
            for row, (mask, ca, cb) in zip(rows, expected):
                selected = [ai[j] for j in ca] + [bi[j] for j in cb]
                assert row["mode"] == "mixed" and (row["u"], row["v"]) == (u, v)
                assert row["selected_local_mask_decimal"] == str(mask)
                assert row["selected_group_indices"] == selected
                assert [entry["group_index"] for entry in row["selected_local_pairs"]] == selected
                assert len(row["Ablock"]["hi"]) == len(row["Ablock"]["lo"]) == 2 * (u + v) + 1
                assert len(row["Bblock"]["hi"]) == len(row["Bblock"]["lo"]) == 2 * (a + b - u - v) + 1
        terminal_rows = cached_block["options"][str(terminal)]
        assert [row["mode"] for row in terminal_rows] == ["all_a", "all_b"]
        assert terminal_rows[0]["selected_group_indices"] == ai
        assert terminal_rows[1]["selected_group_indices"] == bi
        layouts.setdefault(terminal, []).append(2)
    # Canonical enumeration fixes first terminal block to all_a, but no other E.
    total = sum(int(np.prod(layouts[E])) for E in sorted(balance)) + 4
    expected_total = franel(a) * franel(b) // 2 + 2
    assert total == expected_total == raw["balanced_enumeration"]["nonconstant_complement_representatives"]
    # Cached all_a blocks independently reproduce H(x^a) and H(x^b).
    A, B = np.array([1.0]), np.array([1.0])
    for block in cache["blocks"]:
        fair = block["options"][str(terminal)][0]
        A = np.convolve(A, unpair(fair["Ablock"]))
        B = np.convolve(B, unpair(fair["Bblock"]))
    A /= A[0]
    B /= B[0]
    exact_A = np.asarray(substitute(h, a), dtype=np.float64)
    exact_B = np.asarray(substitute(h, b), dtype=np.float64)
    fair_error = max(float(np.max(np.abs(A - exact_A))), float(np.max(np.abs(B - exact_B))))
    assert fair_error < 1e-9
    report = {"status": "passed_numerical_input_and_layout_checks_only", "candidate": raw["name"], "a": a, "b": b,
              "integer_binary_C": True, "block_count": 3, "total_groups": len(raw["groups"]),
              "E_options": expected_E, "radices_by_E": {str(E): layouts[E] if E != terminal else [1, 2, 2] for E in expected_E},
              "canonical_total": total, "franel_formula_total": expected_total, "cached_fair_all_a_max_abs_error": fair_error,
              "input_sha256": sha(args.input), "options_sha256": sha(args.options), "verifier_sha256": sha(Path(__file__)),
              "claim_boundary": "Numerical root/cache validation only; no exact sign exclusion or counterexample is certified."}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
