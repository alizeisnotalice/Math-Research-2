#!/usr/bin/env python3
"""Compare normalized prefix ledgers and independently replay the root and tree checks."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import gcd
from pathlib import Path

from s2_independent_prefix_replay import (
    O, Z, add, atomic, ci, conj, convcoef, cores, neg, sc, shifts, sign, unit,
)

POINTS = ((1, 1), (1, 2), (2, 1), (1, 4), (4, 1),
          (3, 4), (4, 3), (6, 5), (5, 6))
SAMPLE_POWERS = [
    tuple(p**k * q**(16-k) for k in range(17)) for p, q in POINTS
]


def load(path):
    return json.loads(path.read_text())


def core_files(folder):
    return {p.name for p in folder.glob("L*_T*_N*.json")}


def encoded(row):
    return json.dumps(row, separators=(",", ":")).encode() + b"\n"


def row4_production(row):
    L, T, N, D, b, c, g1, q1, rho, sigma, eta, zeta = row
    return L, T, N, D, b, c, g1, q1, rho-L*g1, sigma-L*q1, eta-b, zeta-b


def row58_production(row):
    (L, T, N, D, b, c, d, e, f, g1, q1, rho, sigma, eta, zeta, eta5,
     zeta5, g6, q6, eta7, zeta7) = row
    return (L, T, N, D, b, c, d, e, f, g1, q1, rho-L*g1, sigma-L*q1,
            eta-b, zeta-b, eta5-d, zeta5-d, g6, q6, eta7-f, zeta7-f)


def norm_polynomial(row):
    L, T, N, D, b, c, d, e, f = row
    twice_f_plus_t = (2, 2*L, T, 2*b, 2*c, 2*d, 2*e, 2*f, 2)
    result = [0] * 17
    for i, coeff in enumerate(twice_f_plus_t):
        for j, other in enumerate(twice_f_plus_t):
            result[i+j] += coeff * other
    result[4] -= D
    assert result[0] == result[-1] == 4
    assert D == T*T - 4*N and all(value % 4 == 0 for value in result)
    return tuple(value // 4 for value in result)


def schur_ascending(coeffs):
    """Root count in the open unit disk via real reciprocal subtraction."""
    n = len(coeffs) - 1
    while n:
        c, a = coeffs[0], coeffs[n]
        if abs(c) >= abs(a):
            return n
        coefficients = [a*coeffs[i+1] - c*coeffs[n-i-1] for i in range(n)]
        content = gcd(*coefficients)
        assert content > 0
        coeffs = [x // content for x in coefficients]
        n -= 1
    return None


def root_screen(row):
    h = norm_polynomial(row)
    for j, weights in enumerate(SAMPLE_POWERS):
        if sum(v * w for v, w in zip(h, weights)) <= 0:
            return (10 + j, h)
    outer = schur_ascending([value << k for k, value in enumerate(h)])
    if outer is not None:
        return (30 + outer, h)
    inner = schur_ascending([value << k for k, value in enumerate(reversed(h))])
    if inner is not None:
        return (50 + inner, h)
    return (255, h)


def one_core(job):
    name, expected_core, local, ref, old = job
    result = {}
    for stage, trans in (("stage4", row4_production), ("stage58", row58_production)):
        local_data = load(local/stage/name)
        ref_data = load(ref/stage/"cores"/name)
        assert local_data["status"] == ref_data["status"] == "complete"
        assert tuple(local_data["core"]) == tuple(ref_data["binding"]["core"]) == expected_core
        independent_rows = Counter(map(tuple, local_data["records"]))
        production_rows = Counter(trans(row) for row in ref_data["records"])
        assert independent_rows == production_rows, (stage, name, independent_rows-production_rows,
                                                         production_rows-independent_rows)
        result[stage] = len(local_data["records"])
        if stage == "stage4":
            old_data = load(old/"prefix-full-rows1to4-v1"/"cores"/name)
            old_rows = Counter(trans(row) for row in old_data["records"])
            result["old_extra"] = sorted((old_rows-independent_rows).elements())
            result["old_missing"] = sorted((independent_rows-old_rows).elements())
            result["old_common"] = sum((old_rows & independent_rows).values())
            reference_counters = ref_data["summary"]["counts"]
            for k in ("row2", "row3", "row4"):
                assert local_data["counts"].get(k, 0) == reference_counters.get(k+"_states", 0)
            result["rows4"] = sorted(independent_rows.elements())
        else:
            result["rows58"] = sorted(independent_rows.elements())
    expected_rows = {row[:9] for row in result["rows58"]}
    root_data = load(ref/"root"/"cores"/name)
    assert root_data["status"] == "complete" and tuple(root_data["binding"]["core"]) == expected_core
    reference = {tuple(record["row"]): (record["code"], tuple(record["H"]))
                 for record in root_data["records"]}
    assert expected_rows == reference.keys() and len(root_data["records"]) == len(reference)
    codes = Counter()
    passed = []
    for row in sorted(expected_rows):
        trial = root_screen(row)
        assert trial == reference[row], (name, row, trial, reference[row])
        codes[trial[0]] += 1
        if trial[0] == 255:
            passed.append(row)
    result.update(codes=codes, passed=passed, mothers=len(expected_rows))
    return result


def base_product(poly, prefix, k):
    base = Z
    for j in range(1, min(len(poly)-1, k)+1):
        base = add(base, sc(prefix[k-j], poly[j]))
    return base


def allowed_k(base, k, D):
    if k == 4:
        return {0: base} if sign(base, D) >= 0 and sign(add(O, neg(base)), D) >= 0 else {}
    if k == 6:
        return {-base[0]//2: Z} if base[1] == 0 and base[0] % 2 == 0 else {}
    return dict(shifts(base, D))


def tree(row, depth=14):
    L, T, N, D, b, c, d, e, f = row
    P = [O, ci(L), (T, 1), ci(b), ci(c), ci(d), ci(e), ci(f), O]
    Pb = list(map(conj, P))
    H = norm_polynomial(row)
    layer = {((1,), (1,)): ((1,), (O,), (O,))}
    trace = []
    for k in range(1, depth+1):
        children = {}
        rejected = 0
        for (G, Q), (W, A, B) in layer.items():
            left = allowed_k(base_product(P, G, k), k, D)
            right = allowed_k(base_product(Pb, Q, k), k, D)
            cross_w = sum(G[i] * Q[k-i] for i in range(1, k))
            cross_h = sum(H[j] * W[k-j] for j in range(1, min(16, k)+1))
            accepted = 0
            for gk, ak in left.items():
                for ck in (0, 1):
                    wk = ck - cross_h
                    qk = wk - gk - cross_w
                    if qk not in right:
                        continue
                    bk = right[qk]
                    AA, BB = A + (ak,), B + (bk,)
                    assert convcoef(AA, BB, k, D) == ci(ck)
                    child = (G+(gk,), Q+(qk,))
                    assert child not in children
                    children[child] = (W+(wk,), AA, BB)
                    accepted += 1
            rejected += len(left) * len(right) - accepted
        statehash = hashlib.sha256()
        for child in sorted(children):
            statehash.update(encoded(child))
        trace.append({"k": k, "live": len(children), "nonbinary_C": rejected,
                      "state_sha256": statehash.hexdigest()})
        if not children:
            return {"status": "dead_all_prefixes", "depth_reached": k, "trace": trace}
        layer = children
    return {"status": "unknown_depth", "depth_reached": depth, "trace": trace}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--local", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--old", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--trees", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    expected = {f"L{L}_T{T}_N{N}.json": (L,T,N,D) for L,T,N,D in cores()}
    assert len(expected) == 11760
    for stage in ("stage4", "stage58"):
        assert core_files(args.local/stage) == core_files(args.reference/stage/"cores") == expected.keys()
    assert core_files(args.reference/"root"/"cores") == expected.keys()
    assert core_files(args.old/"prefix-full-rows1to4-v1"/"cores") == expected.keys()
    jobs = [(name, core, args.local, args.reference, args.old) for name, core in sorted(expected.items())]
    rows4, rows58, passed, extra, missing = [], [], [], [], []
    codes, counts = Counter(), Counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(one_core, jobs, chunksize=16):
            codes.update(result["codes"])
            for key in ("stage4", "stage58", "old_common", "mothers"):
                counts[key] += result[key]
            rows4.extend(result["rows4"])
            rows58.extend(result["rows58"])
            passed.extend(result["passed"])
            extra.extend(result["old_extra"])
            missing.extend(result["old_missing"])
    passed.sort()
    assert len(passed) == len(set(passed)) == codes[255]
    counts.update(old_extra=len(extra), old_missing=len(missing), root_pass=len(passed))
    atomic(args.out/"prefix-differences.json", {"old_extra": sorted(extra), "old_missing": sorted(missing)})
    atomic(args.out/"passed.json", passed)
    trace_hist, max_live_hist = Counter(), Counter()
    deepest = []
    trace_layers = 0
    if args.trees is not None:
        assert len(list(args.trees.glob("mother-*.json"))) == len(passed)
    for row in passed:
        result = tree(row, depth=30)
        assert result["status"] == "dead_all_prefixes", (row, result)
        trace_hist[result["depth_reached"]] += 1
        trace_layers += result["depth_reached"]
        max_live_hist[max(layer["live"] for layer in result["trace"])] += 1
        identifier = hashlib.sha256(encoded(row).rstrip()).hexdigest()
        path = f"mother-{identifier}.json"
        if args.trees is not None:
            ref = load(args.trees/path)
            assert tuple(ref["binding"]["row"]) == row and ref["status"] == "complete"
            assert result == ref["ab"] == ref["hgq"], (row, result, ref)
        if not deepest or result["depth_reached"] > deepest[0][1]:
            deepest = [(row, result["depth_reached"])]
        elif result["depth_reached"] == deepest[0][1]:
            deepest.append((row, result["depth_reached"]))
        atomic(args.out/"trees"/path, {"row": row, **result})
    summary = dict(status="complete", cores=len(expected), counts=dict(counts), root_codes=dict(codes),
                   tree_deaths=dict(trace_hist), tree_layers=trace_layers, max_live_hist=dict(max_live_hist),
                   deepest=deepest, production_trees_compared=args.trees is not None,
                   sha256={name: hashlib.sha256(b"".join(map(encoded, sorted(rows)))).hexdigest()
                           for name, rows in (("row4", rows4), ("row58", rows58),
                                              ("mothers", [row[:9] for row in rows58]),
                                              ("root_pass", passed))})
    atomic(args.out/"summary.json", summary)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
