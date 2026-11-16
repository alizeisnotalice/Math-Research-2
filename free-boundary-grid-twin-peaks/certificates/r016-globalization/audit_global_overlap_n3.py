#!/usr/bin/env python3
"""Exact n=3 audit of cross-state overlap in a canonical gate construction.

This is only a finite audit of the multiplicity identity in
globalization_accounting.md.  It takes every canonical ordered partition
V=A dotcup B dotcup C with connected non-adjacent nonempty A,B,C, assigns
each C vertex its lexicographically first A-union-B neighbor, selects all C,
and emits the bounded-parent gate family.  It then counts how many different
partitions emit each resulting two-birth word.
"""

from collections import Counter
from itertools import combinations, permutations
from math import ceil, comb, factorial


N = 3
V = tuple((x, y) for x in range(N) for y in range(N))
INDEX = {v: i for i, v in enumerate(V)}


def adj(v):
    x, y = v
    return {(x + dx, y + dy) for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))
            if 0 <= x + dx < N and 0 <= y + dy < N}


ADJ = {v: adj(v) for v in V}


def connected(S):
    if not S:
        return False
    seen = {next(iter(S))}
    todo = list(seen)
    while todo:
        v = todo.pop()
        for u in ADJ[v] & S:
            if u not in seen:
                seen.add(u)
                todo.append(u)
    return seen == S


def onebirth_words(S):
    ans = []
    for word in permutations(sorted(S)):
        seen = set()
        for v in word:
            if seen and not (ADJ[v] & seen):
                break
            seen.add(v)
        else:
            ans.append(word)
    return ans


def root_first_interleavings(x, y):
    # Both first letters occur in global positions 0,1; internal orders fixed.
    for first, second in ((x, y), (y, x)):
        for slots in combinations(range(len(x) + len(y) - 2), len(x) - 1):
            slots = set(slots)
            ix = iy = 1
            out = [first[0], second[0]]
            for slot in range(len(x) + len(y) - 2):
                if slot in slots:
                    out.append(x[ix]); ix += 1
                else:
                    out.append(y[iy]); iy += 1
            yield tuple(out)


def gate_extensions(alpha, children, parents):
    """All gate-tree extensions for the degree-four capacity construction."""
    M, k = len(alpha), len(children)
    q = ceil(k / 4)
    r = k - 4 * (q - 1)
    assert q <= M - 1
    rank = {v: i + 1 for i, v in enumerate(alpha)}
    children = tuple(sorted(children, key=lambda v: (rank[parents[v]], v)))
    for i, child in enumerate(children, 1):
        assert rank[parents[child]] <= M - ((k - i) // 4)
    prefix = alpha[:M - q]
    gates = alpha[M - q:]
    groups = [children[:r]]
    cursor = r
    for _ in range(1, q):
        groups.append(children[cursor:cursor + 4])
        cursor += 4
    nodes = tuple(gates) + children
    pred = {v: set() for v in nodes}
    for j, gate in enumerate(gates):
        if j:
            pred[gate].add(gates[j - 1])
        for child in groups[j]:
            pred[child].add(gate)
    ans = []
    def rec(word, seen):
        if len(word) == len(nodes):
            ans.append(tuple(word)); return
        for v in nodes:
            if v not in seen and pred[v] <= seen:
                rec(word + [v], seen | {v})
    rec([], set())
    expected = factorial(k + q - 1) // (5 ** (q - 1) * factorial(q - 1))
    assert len(ans) == expected
    return prefix, ans


def births(word):
    seen, total = set(), 0
    for v in word:
        if not (ADJ[v] & seen):
            total += 1
        seen.add(v)
    return total


def all_partitions():
    # Canonical orientation: min(A)<min(B), so A/B swapping is not counted.
    for colors in range(3 ** len(V)):
        z = colors
        A, B, C = set(), set(), set()
        for v in V:
            digit, z = z % 3, z // 3
            (A if digit == 0 else B if digit == 1 else C).add(v)
        if not A or not B or not C or min(A) > min(B):
            continue
        if not connected(A) or not connected(B):
            continue
        if any(ADJ[v] & B for v in A):
            continue
        if any(not (ADJ[v] & (A | B)) for v in C):
            continue
        M, k = len(A) + len(B), len(C)
        if ceil(k / 4) > M - 1:
            continue
        parents = {v: min(ADJ[v] & (A | B)) for v in C}
        yield A, B, C, parents


def main():
    multiplicity = Counter()
    state_count = output_sum = 0
    max_within = 0
    for A, B, C, parents in all_partitions():
        emitted = set()
        for aw in onebirth_words(A):
            for bw in onebirth_words(B):
                for alpha in root_first_interleavings(aw, bw):
                    prefix, tails = gate_extensions(alpha, tuple(C), parents)
                    for tail in tails:
                        word = prefix + tail
                        assert births(word) == 2
                        emitted.add(word)
        # The fixed-state recovery proof predicts this cardinality exactly.
        a, b, k = len(A), len(B), len(C)
        q = ceil(k / 4)
        expected = (len(onebirth_words(A)) * len(onebirth_words(B))
                    * 2 * comb(a + b - 2, a - 1)
                    * factorial(k + q - 1) // (5 ** (q - 1) * factorial(q - 1)))
        assert len(emitted) == expected
        state_count += 1
        output_sum += len(emitted)
        max_within = max(max_within, len(emitted))
        for word in emitted:
            multiplicity[word] += 1
    assert output_sum == sum(multiplicity.values())
    assert all(births(word) == 2 for word in multiplicity)
    hist = Counter(multiplicity.values())
    print({
        'n': N,
        'canonical_partitions': state_count,
        'sum_fixed_state_outputs': output_sum,
        'distinct_two_birth_outputs': len(multiplicity),
        'max_cross_state_multiplicity': max(multiplicity.values(), default=0),
        'histogram_multiplicity': dict(sorted(hist.items())),
        'max_single_state_output_size': max_within,
        'double_count_identity_verified': True,
    })


if __name__ == '__main__':
    main()
