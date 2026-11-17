#!/usr/bin/env python3
"""Exact RM(2,4)/G extension classification for the D6 residual.

G|_R is RM(1,4)+<v>, so the quotient of RM(2,4) has 2^5 classes.
For each class h we test f0+h+g bent for every g in G.  A linear quotient
subspace contained in the good classes gives a simultaneous fixed-atom
quadratic residual extension.
"""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path


EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
F0 = (1 << 2) | (1 << 3)       # x0x3+x1x2
VQ = (1 << 1) | (1 << 3) | (1 << 4)  # x0x2+x1x2+x1x3
FREE = [0, 1, 2, 3, 5]         # complement to <VQ>, whose pivot is bit 4


def fwht(a):
    a = a[:]
    h = 1
    while h < len(a):
        for p in range(0, len(a), 2 * h):
            for i in range(p, p + h):
                l, r = a[i], a[i + h]
                a[i], a[i + h] = l + r, l - r
        h *= 2
    return a


def qval(mask, x):
    return sum(((mask >> j) & 1) * (((x >> a) & 1) * ((x >> b) & 1))
               for j, (a, b) in enumerate(EDGES)) & 1


def aval(mask, x):
    # mask bit 0 is constant, bits 1..4 are x0,...,x3 coefficients.
    return ((mask & 1) ^ (((mask >> 1) & x).bit_count() & 1)) & 1


def bent_quad(mask):
    W = fwht([1 if qval(mask, x) == 0 else -1 for x in range(16)])
    return {abs(z) for z in W} == {4}


def lift(qclass):
    return sum(((qclass >> i) & 1) << bit for i, bit in enumerate(FREE))


def qcoord(mask):
    assert ((mask >> 4) & 1) == 0
    return sum(((mask >> bit) & 1) << i for i, bit in enumerate(FREE))


def all_subspaces_contained(good):
    spaces = {frozenset({0})}
    frontier = [frozenset({0})]
    universe = range(32)
    while frontier:
        S = frontier.pop()
        for x in universe:
            if x in S:
                continue
            T = frozenset(y ^ (x if bit else 0) for y in S for bit in (0, 1))
            if len(T) == 2 * len(S) and T <= good and T not in spaces:
                spaces.add(T)
                frontier.append(T)
    d = max(len(S).bit_length() - 1 for S in spaces)
    return d, [S for S in spaces if len(S) == 1 << d], spaces


def h_weight(mask):
    return sum(qval(mask, x) for x in range(16))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, help="JSON result path")
    args = parser.parse_args()
    # all 64 elements of RM(1,4)+<v> represented by affine mask and c.
    G = [(aff, c) for aff in range(32) for c in range(2)]
    good = set()
    detail = []
    for cls in range(32):
        h = lift(cls)
        checks = []
        for aff, c in G:
            vals = [f0bit ^ qval(h, x) ^ aval(aff, x) ^ (c * qval(VQ, x))
                    for x in range(16)
                    for f0bit in [qval(F0, x)]]
            W = fwht([1 if z == 0 else -1 for z in vals])
            checks.append({abs(z) for z in W} == {4})
        if all(checks):
            good.add(cls)
        residual_weights = []
        for aff, c in G:
            residual_weights.append(2 * sum(
                qval(h, x) ^ aval(aff, x) ^ (c * qval(VQ, x))
                for x in range(16)
            ))
        histogram = {str(weight): residual_weights.count(weight)
                     for weight in sorted(set(residual_weights))}
        detail.append({
            "quotient_coordinate": cls,
            "quadratic_edge_mask": h,
            "quadratic_edges": [list(EDGES[j]) for j in range(6) if (h >> j) & 1],
            "all_64_G_translates_bent": all(checks),
            "representative_h_weight_on_F2_4": h_weight(h),
            "residual_word_weight_on_R": 2 * h_weight(h),
            "h_plus_G_residual_weight_histogram": histogram,
            "h_plus_G_contains_q_weight": 8 in residual_weights,
        })
    dmax, maximal, spaces = all_subspaces_contained(good)
    assert len(good) == 6
    assert dmax == 1
    assert len(maximal) == 5
    assert all(0 in S for S in maximal)
    assert all(len(S) == 2 for S in maximal)
    assert all(x["h_plus_G_contains_q_weight"] for x in detail
               if x["all_64_G_translates_bent"] and x["quotient_coordinate"] != 0)
    out = {
        "method": "Exact integer FWHT enumeration over RM(2,4)/(RM(1,4)+<v>).",
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "f0_edge_mask": F0, "v_edge_mask": VQ, "edge_order": EDGES,
        "quotient_dimension": 5,
        "good_class_count": len(good),
        "good_classes": [x for x in detail if x["all_64_G_translates_bent"]],
        "maximum_simultaneous_extension_dimension": dmax,
        "maximum_subspace_count": len(maximal),
    }
    # More readable explicit maximal subspaces.
    out["maximum_subspaces"] = [[
        {
            "quotient_coordinate": x,
            "quadratic_edge_mask": lift(x),
            "quadratic_edges": [list(EDGES[j]) for j in range(6) if (lift(x) >> j) & 1],
            "residual_word_weight_on_R": 2 * h_weight(lift(x)),
        } for x in sorted(S)
    ] for S in maximal]
    out["all_good_detail"] = detail
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print(json.dumps({
        "good_class_count": len(good),
        "maximum_simultaneous_extension_dimension": dmax,
        "maximum_subspace_count": len(maximal),
        "nonzero_good_coset_weight_histograms": {
            str(x["quadratic_edge_mask"]): x["h_plus_G_residual_weight_histogram"]
            for x in detail if x["all_64_G_translates_bent"] and x["quotient_coordinate"]
        },
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
