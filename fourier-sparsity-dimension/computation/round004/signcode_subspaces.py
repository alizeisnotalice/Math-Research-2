#!/usr/bin/env python3
"""Exact sign-code/subspace audit for all Boolean functions on F_2^n, n<=4."""
import argparse
import hashlib
import itertools
import json
import platform
import sys
import time

def fwht(values):
    values = values[:]
    step = 1
    while step < len(values):
        for start in range(0, len(values), 2 * step):
            for j in range(start, start + step):
                left, right = values[j], values[j + step]
                values[j], values[j + step] = left + right, left - right
        step *= 2
    return values

def rref(vectors):
    pivots = {}
    for vector in vectors:
        x = vector
        while x:
            bit = x.bit_length() - 1
            if bit in pivots:
                x ^= pivots[bit]
            else:
                pivots[bit] = x
                break
    return pivots

def reduce_by(vector, pivots):
    x = vector
    for bit in sorted(pivots, reverse=True):
        if (x >> bit) & 1:
            x ^= pivots[bit]
    return x

def span(basis):
    values = {0}
    for vector in basis:
        values |= {x ^ vector for x in values}
    return values

def compact_map(c_pivots, support_size):
    free = [bit for bit in range(support_size) if bit not in c_pivots]
    def compact(vector):
        remainder = reduce_by(vector, c_pivots)
        return sum(((remainder >> bit) & 1) << index for index, bit in enumerate(free))
    return free, compact

def max_subspace(q_members, quotient_dimension):
    """Exhaust all linear W contained in q_members by closing candidate spans."""
    q_members = set(q_members)
    if 0 not in q_members:
        raise AssertionError("quotient must contain zero")
    if len(q_members) == 1 << quotient_dimension:
        return quotient_dimension, [1 << i for i in range(quotient_dimension)], 1
    seen = {frozenset({0})}
    best_basis = []
    examined = 0
    def dfs(basis, current):
        nonlocal best_basis, examined
        examined += 1
        if len(basis) > len(best_basis):
            best_basis = basis[:]
        for vector in q_members:
            if vector == 0 or vector in current:
                continue
            enlarged = current | {x ^ vector for x in current}
            if not enlarged <= q_members:
                continue
            key = frozenset(enlarged)
            if key in seen:
                continue
            seen.add(key)
            dfs(basis + [vector], enlarged)
    dfs([], {0})
    return len(best_basis), best_basis, examined

def verify_graph(coefficients, support, lifts, n):
    """Check that every z-slice of the lifted graph has a ±1 inverse Walsh table."""
    r = len(lifts)
    raw_scale = 1 << n
    for z in range(1 << r):
        transformed = [0] * (1 << n)
        for index, frequency in enumerate(support):
            sign = -1 if sum(((z >> j) & 1) * ((lifts[j] >> index) & 1) for j in range(r)) % 2 else 1
            transformed[frequency] = sign * coefficients[frequency]
        truth_raw = fwht(transformed)
        if any(abs(value) != raw_scale for value in truth_raw):
            return False
    return True

def inspect_group(n, group):
    _, coefficients = group[0]
    support = [frequency for frequency, value in enumerate(coefficients) if value]
    size = len(support)
    relative = []
    for _, other in group:
        mask = 0
        for index, frequency in enumerate(support):
            if other[frequency] != coefficients[frequency]:
                mask |= 1 << index
        relative.append(mask)
    k = set(relative)
    c_generators = [(1 << size) - 1]
    for coordinate in range(n):
        c_generators.append(sum(((frequency >> coordinate) & 1) << index for index, frequency in enumerate(support)))
    c_pivots = rref(c_generators)
    c_basis = [c_pivots[bit] for bit in sorted(c_pivots, reverse=True)]
    c_all = span(c_basis)
    if not c_all <= k:
        raise AssertionError("C is not contained in K")
    if any((code ^ translate) not in k for code in k for translate in c_all):
        raise AssertionError("K is not C-invariant")
    free, compact = compact_map(c_pivots, size)
    quotient_to_rep = {}
    for code in k:
        quotient_to_rep.setdefault(compact(code), code)
    q_members = set(quotient_to_rep)
    if len(k) != len(q_members) * len(c_all):
        raise AssertionError("incomplete C coset")
    r, q_basis, examined = max_subspace(q_members, len(free))
    lifts = [quotient_to_rep[vector] for vector in q_basis]
    if not verify_graph(coefficients, support, lifts, n):
        raise AssertionError("lifted graph is not Boolean")
    return {"n": n, "s": size, "group_size": len(group), "dimC": len(c_basis),
            "quotient_dimension": len(free), "quotient_members": len(q_members),
            "subspaces_examined": examined, "r": r, "dimV": len(c_basis) + r,
            "support": support, "integercoeffs": [coefficients[frequency] for frequency in support],
            "C_basis": c_basis, "lifted_signcode_basis": lifts,
            "graph_boolean_verified": True}

def enumerate_n(n):
    points = 1 << n
    groups = {}
    for truth_mask in range(1 << points):
        values = [-1 if (truth_mask >> x) & 1 else 1 for x in range(points)]
        coefficients = tuple(fwht(values))
        groups.setdefault(tuple(abs(value) for value in coefficients), []).append((truth_mask, coefficients))
    return groups

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=4)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    began = time.monotonic()
    report = {"method": "complete_integer_truth_table_walsh_signcode_quotient_subspaces", "max_n": args.max_n,
              "by_n": {}, "by_s": {}}
    for n in range(args.max_n + 1):
        groups = enumerate_n(n)
        rows = []
        for group in groups.values():
            row = inspect_group(n, group)
            rows.append(row)
            old = report["by_s"].get(str(row["s"]))
            if old is None or row["dimV"] > old["dimV"]:
                report["by_s"][str(row["s"])] = row
        report["by_n"][str(n)] = {"functions": 1 << (1 << n), "groups": len(groups), "rows": rows,
                                    "max_quotient_dimension": max(row["quotient_dimension"] for row in rows)}
    report["elapsed_seconds"] = time.monotonic() - began
    report["python"] = sys.version
    report["platform"] = platform.platform()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    open(args.out, "w").write(text)
    print(json.dumps({"status": "passed", "sha256": hashlib.sha256(text.encode()).hexdigest(),
                      "elapsed_seconds": report["elapsed_seconds"],
                      "max_quotient_dimension": max(x["max_quotient_dimension"] for x in report["by_n"].values())}))

if __name__ == "__main__":
    main()
