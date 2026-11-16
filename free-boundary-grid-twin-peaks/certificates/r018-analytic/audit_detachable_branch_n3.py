#!/usr/bin/env python3
"""Exact n=3 audit of b=1 detachable-branch insertions and collisions."""

from collections import Counter
from itertools import combinations, permutations
from math import comb
import json
from pathlib import Path


def graph(n):
    adj = []
    for y in range(n):
        for x in range(n):
            ns = []
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    ns.append(yy * n + xx)
            adj.append(tuple(ns))
    return adj


def peaks(order, adj):
    pos = [0] * len(order)
    for i, v in enumerate(order):
        pos[v] = i
    return tuple(v for v, ns in enumerate(adj) if all(pos[v] < pos[u] for u in ns))


def is_onepeak(order, adj):
    return len(peaks(order, adj)) == 1


def shuffles(left, right):
    for slots in combinations(range(len(left) + len(right)), len(left)):
        slots = set(slots)
        i = j = 0
        out = []
        for position in range(len(left) + len(right)):
            if position in slots:
                out.append(left[i])
                i += 1
            else:
                out.append(right[j])
                j += 1
        yield tuple(out)


def main():
    adj = graph(3)
    sources = 0
    events = 0
    proposed_outputs = 0
    output_load = Counter()
    invalid_outputs = 0
    canonical_weight_histogram = Counter()
    zero_example = None
    # L=A,w,v,T; A is the first a vertices and v has no neighbour in A.
    for order in permutations(range(9)):
        if not is_onepeak(order, adj):
            continue
        sources += 1
        branch_weight = 0
        for a in range(1, 7):
            A = order[:a]
            w, v = order[a], order[a + 1]
            if v in adj[w] and all(v not in adj[x] for x in A):
                events += 1
                for prefix in shuffles(A, (v,)):
                    output = prefix + (w,) + order[a + 2:]
                    proposed_outputs += 1
                    if len(peaks(output, adj)) != 2:
                        invalid_outputs += 1
                    output_load[output] += 1
            # General detachable branch B immediately after w in the source word.
            for b in range(1, 9 - a):
                B = order[a + 1:a + 1 + b]
                if min(A) >= min(B):  # canonical orientation for unique recovery
                    continue
                if any(y in adj[x] for x in A for y in B):
                    continue
                if B[0] not in adj[w]:
                    continue
                if any(not any(B[i] in adj[B[j]] for j in range(i)) for i in range(1, b)):
                    continue
                branch_weight += comb(a + b, a)
        canonical_weight_histogram[branch_weight] += 1
        if branch_weight == 0 and zero_example is None:
            zero_example = list(order)
    ans = {
        "A1_n3": sources,
        "source_events_b1": events,
        "proposed_outputs": proposed_outputs,
        "invalid_outputs": invalid_outputs,
        "distinct_outputs": len(output_load),
        "max_output_load": max(output_load.values(), default=0),
        "load_histogram": dict(sorted(Counter(output_load.values()).items())),
        "canonical_detachable_branch_weight_histogram": dict(sorted(canonical_weight_histogram.items())),
        "canonical_detachable_branch_total_weight": sum(
            weight * multiplicity for weight, multiplicity in canonical_weight_histogram.items()
        ),
        "canonical_detachable_branch_zero_weight_example": zero_example,
        "status": "passed" if invalid_outputs == 0 else "failed",
    }
    Path("detachable-branch-n3.json").write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
