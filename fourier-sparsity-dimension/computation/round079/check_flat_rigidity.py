"""Run with python3 -B check_flat_rigidity.py. All arithmetic is integral."""
import itertools
import json
import platform
from collections import Counter
from functools import lru_cache


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


def rank(vectors):
    pivots = {}
    for x in vectors:
        while x and x.bit_length() in pivots:
            x ^= pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()] = x
    return len(pivots)


@lru_cache(None)
def components(labels):
    columns = [1+2*v for v in labels]
    r = rank(columns)
    parts = [{i} for i in range(len(columns))]
    for size in range(2, len(columns)+1):
        for choice in itertools.combinations(range(len(columns)), size):
            subset = [columns[i] for i in choice]
            if rank(subset) == size or any(rank(subset[:i]+subset[i+1:]) < size-1
                                           for i in range(size)):
                continue
            involved = [part for part in parts if set(choice) & part]
            parts = [part for part in parts if not set(choice) & part] + [set.union(*involved)]
    return r, [[labels[i] for i in sorted(part)] for part in parts]


def all_flats(n):
    spaces = {frozenset({0})}
    for dim in range(n+1):
        yield n-dim, {tuple(sorted(x ^ b for x in h)) for h in spaces for b in range(2**n)}
        spaces = {h | frozenset(x ^ b for x in h)
                  for h in spaces for b in range(2**n) if b not in h}


def audit_tables(n):
    N = 2**n
    flats = list(all_flats(n))
    count = Counter()
    equal = Counter()
    types = Counter()
    for bits in range(2**N):
        transformed = walsh([1-2*((bits >> x) & 1) for x in range(N)])
        for t, planes in flats:
            K = 2**(t+1)
            for plane in planes:
                labels = tuple(e for e in plane if transformed[e])
                r, parts = components(labels)
                q = len(labels)
                nullity = q-r
                p = sum(len(part) == 1 for part in parts)
                assert q+p <= 2*nullity+K
                form = True
                circuit_capacity = 0
                for part in parts:
                    coeffs = [transformed[e] for e in part]
                    if len(part) == 1:
                        form &= K*abs(coeffs[0]) == 2*N
                    else:
                        form &= len(part) >= 4 and len(part) % 2 == 0
                        form &= rank([1+2*e for e in part]) == len(part)-1
                        form &= all(K*abs(a) == N for a in coeffs)
                        form &= sum(a < 0 for a in coeffs) % 2 == 1
                        circuit_capacity += len(part)-2
                form &= circuit_capacity+2*p == K
                assert bool(form) == (q+p == 2*nullity+K)
                count[t] += 1
                if form:
                    equal[t] += 1
                    types[f"{t},{p},{nullity}"] += 1
    return {"dimension": n, "table_flat_pairs_by_t": dict(count),
            "equality_by_t": dict(equal), "equality_by_t_p_nullity": dict(types)}


def check_constructions():
    totals = {}
    for t in range(5):
        K = 2**(t+1)
        count = 0
        for p in range(K//2+1):
            for l in range(K//2-p+1):
                if not l and 2*p != K:
                    continue
                lengths = [4]*l
                if l:
                    lengths[-1] += K-2*p-2*l
                sums = {2*(2*i-p) for i in range(p+1)}
                for length in lengths:
                    sums = {x+length-2*j for x in sums for j in range(1, length, 2)}
                assert sums == set(range(-K, K+1, 4))
                for x in sums:
                    positives = (K+x)//4
                    assert 2*positives-2**t == x//2
                assert sum(lengths)+2*p == 2*l+K
                count += 1
        totals[t] = count
    return totals


if __name__ == "__main__":
    print(json.dumps({"python": platform.python_version(), "seed": None, "tolerance": 0,
                      "three_variable_exhaustion": audit_tables(3),
                      "sharp_construction_parameter_cases_by_t": check_constructions()}, indent=2))
