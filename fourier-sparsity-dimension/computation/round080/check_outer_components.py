"""Run python3 -B check_outer_components.py. Exact checks on all three-bit tables."""
from collections import Counter, defaultdict
from functools import lru_cache
import json
import platform


def parity(x):
    return x.bit_count() & 1


def walsh(values):
    out = list(values)
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2*width):
            for j in range(start, start+width):
                a, b = out[j], out[j+width]
                out[j], out[j+width] = a+b, a-b
        width *= 2
    return out


def extend(subspace, v):
    return frozenset(subspace | {x ^ v for x in subspace})


def span(labels):
    out = frozenset({0})
    for x in labels:
        if x not in out:
            out = extend(out, x)
    return out


@lru_cache(None)
def subspaces(space):
    layer = {frozenset({0})}
    result = []
    while layer:
        result.extend(layer)
        layer = {extend(H, v) for H in layer for v in space-H}
    return result


@lru_cache(None)
def parts(labels):
    # The fundamental circuits of any basis generate the matroid components.
    expression = {0: frozenset()}
    groups = [{v} for v in labels]
    for v in labels:
        if v in expression:
            circ = expression[v] | {v}
            hit = [g for g in groups if g & circ]
            groups = [g for g in groups if not g & circ] + [set.union(*hit)]
        else:
            expression.update({a ^ v: subset | {v} for a, subset in expression.items()})
    return groups


def audit(n):
    N = 2**n
    stats = Counter()
    nested = []
    for table in range(2**N):
        f = [1-2*((table >> x) & 1) for x in range(N)]
        values = f + [-a for a in f]
        coefficients = {N+e:a for e, a in enumerate(walsh(f)) if a}
        E = set(coefficients)
        V = span(E)
        for H in subspaces(V):
            pi = {v:min(v ^ h for h in H) for v in V}
            W = {pi[v] for v in V}
            components = parts(tuple(sorted(E-H)))
            m = len(components)
            stats[f"components_{m}"] += 1
            assert m <= 3
            if m < 2:
                continue
            U = [span(pi[v] for v in C) for C in components]
            if m == 3:
                assert len(W) == 4 and all(len(u) == 2 for u in U)
                for x in range(2*N):
                    bands = Counter()
                    for v, a in coefficients.items():
                        bands[pi[v]] += a*(-1)**parity(v & x)
                    assert all(2*abs(bands[w]) == N for w in W)
                    assert sum(bands[w]<0 for w in W) % 2 == 1
                continue
            if len(U[0]) > len(U[1]):
                U.reverse()
            assert U[0] <= U[1] and len(U[0]) == 2 and U[1] == W
            if len(W) == 2:
                continue
            b = next(iter(U[0]-{0}))
            fibers = defaultdict(lambda: [set(), set()])
            for x, value in enumerate(values):
                key = tuple(parity(h & x) for h in sorted(H))
                fibers[key][parity(b & x)].add(value)
            assert all(min(map(len, pair)) == 1 for pair in fibers.values())
            stats["nested_two_with_rank_at_least_two"] += 1
            if len(nested) < 1:
                nested.append({"truth_table":table, "H":sorted(H), "quotient_size":len(W)})
    return {"n":n, "counters":dict(stats), "nested_example":nested}


def selector(m):
    n = m+2
    N = 2**n
    values = []
    for raw in range(N):
        z, x, y = raw & 1, (raw >> 1) & 1, raw >> 2
        values.append((1-2*x) if z == 0 else (1-2*int(y == 0)))
    E = {N+index for index, value in enumerate(walsh(values)) if value}
    V = span(E)
    H = frozenset({0, N+2})
    assert H <= V
    blocks = parts(tuple(sorted(E-H)))
    U = sorted([span(min(v^h for h in H) for v in block) for block in blocks],key=len)
    assert len(U) == 2 and len(U[0]) == 2 and U[0] < U[1]
    assert len(U[1]) == len(V)//2 == 2**(m+2)
    b = next(iter(U[0]-{0}))
    lifted = values+[-a for a in values]
    fibers = defaultdict(lambda: [set(),set()])
    for x, value in enumerate(lifted):
        fibers[parity((N+2) & x)][parity(b & x)].add(value)
    assert all(sorted(map(len, halves)) == [1,2] for halves in fibers.values())
    return {"data_variables":m, "sparsity":len(E), "quotient_rank":m+2,
            "component_sizes":sorted(map(len,blocks)), "nonconstant_slices_fraction":"1/2"}


if __name__ == "__main__":
    print(json.dumps({"python":platform.python_version(),"seed":None,"tolerance":0,
                      "exhaustion":audit(3), "selectors":[selector(m) for m in range(2,7)]},indent=2))
