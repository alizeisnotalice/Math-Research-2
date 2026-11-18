"""Exact valuation-code checks in F_2; no floating point or dependencies."""
import json
import platform
from collections import Counter
from itertools import product


def fwht(values):
    out = list(values)
    half = 1
    while half < len(out):
        for off in range(0, len(out), 2 * half):
            for u in range(off, off + half):
                a, b = out[u], out[u + half]
                out[u], out[u + half] = a + b, a - b
        half *= 2
    return out


def valuation(x):
    return (x & -x).bit_length() - 1


def span(vectors):
    found = {0}
    for v in vectors:
        if v not in found:
            found |= {u ^ v for u in found}
    return found


def inspect(F):
    n = len(F).bit_length() - 1
    w = fwht(F)
    v0 = min(valuation(z) for z in w if z)
    k = n - v0
    if not k:
        return None
    odd = [a for a, z in enumerate(w) if z and valuation(z) == v0]
    shift = odd[0]
    W = span(x ^ shift for x in odd)
    if len(W) != len(odd):
        return None
    to_coset = {x: min(x ^ u for u in W) for x in range(len(F))}
    tails = {x ^ shift: valuation(z)-v0 for x, z in enumerate(w)
             if z and valuation(z) > v0}
    # Minimal bitmask representatives of linear cosets are a linear section.
    quot = {to_coset[x] for x in tails}
    j_beta = {b: min(tails[x] for x in tails if to_coset[x] == b) for b in quot}
    lower = {0}
    levels = []
    for j in sorted(set(j_beta.values())):
        upper = span(lower | {b for b in quot if j_beta[b] == j})
        dim = (len(upper) // len(lower)).bit_length() - 1
        points = [x for x, level in tails.items() if level == j and to_coset[x] not in lower]
        fibres = Counter(to_coset[x] for x in points)
        assert all(m % 2 == 0 for m in fibres.values())
        weights = set()
        for v in range(len(F)):
            if any((v & x).bit_count() % 2 for x in W | lower):
                continue
            weight = sum((v & x).bit_count() % 2 for x in points)
            if weight:
                weights.add(weight)
        delta = 1 << (k - j - 1)
        bound = 2 * sum((delta + (1 << i) - 1) >> i for i in range(dim))
        assert len(points) >= bound
        if dim:
            assert min(weights) >= 2 * delta
        levels.append({"j": j, "rank_gain": dim, "paired_length": len(points) // 2,
                       "paired_min_distance": min(weights)//2 if weights else None,
                       "point_lower_bound": bound})
        lower = upper
    return {"k": k, "core": len(odd), "tail": len(tails), "levels": levels}


def selectors():
    # t selects an odd-weight flip on y or the independently tagged data phase.
    cases = {}
    for data in ("character", "two_bit_bent"):
        a, d = 4, 2
        F = []
        for index in range(1 << (1+a+d+1)):
            t = index & 1
            y = (index >> 1) & ((1 << a)-1)
            z = (index >> (a+1)) & 3
            tag = index >> (a+d+1)
            phase = tag ^ ((z == 3) if data == "two_bit_bent" else 0)
            F.append((-1)**(phase if t else int(y == 0)))
        result = inspect(F)
        assert sum(e["point_lower_bound"] for e in result["levels"]) == result["tail"]
        cases[data] = result
    return cases


census = {}
for n in range(1, 5):
    count = 0
    with_tail = 0
    for F in product((-1, 1), repeat=1 << n):
        result = inspect(F)
        if result:
            count += 1
            with_tail += bool(result["tail"])
    census[n] = {"flat_core_noncharacter": count, "with_tail": with_tail}
print(json.dumps({"python": platform.python_version(), "complete_census": census,
                  "selectors_attaining_layer_bound": selectors()}, indent=2))
