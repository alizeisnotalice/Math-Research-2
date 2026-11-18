"""Certify the optimal support norm squared for every four-variable Boolean spectrum."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import platform
from support_census import cuts, enlarge


def permute(mask, perm):
    return sum(1 << j for i, j in enumerate(perm) if mask >> i & 1)


def walsh_supports():
    columns = [sum(1 << x for x in range(16) if (x & i).bit_count() & 1)
               for i in range(16)]
    return Counter(sum(1 << i for i, col in enumerate(columns)
                       if (col ^ table).bit_count() != 8) for table in range(32768))


def orbit(seed, generators, action):
    seen = {seed}
    work = [seed]
    while work:
        for image in action(work.pop(), generators):
            if image not in seen:
                seen.add(image)
                work.append(image)
    return seen


def certificate(mask):
    points = [i for i in range(16) if mask >> i & 1]
    span = frozenset([0])
    for i in points:
        span = enlarge(span, i + 16)
    constraints = {(mask & ~inner, demand) for inner, demand in cuts(span)}
    # The optimum is the projection onto its active equalities. Find a positive
    # dual certificate using exact coordinate steps on cut orbits by their size.
    groups = {}
    for outside, demand in constraints:
        groups.setdefault((outside.bit_count(), demand), []).append(outside)
    rows = [[sum((outside >> i) & 1 for outside in group) for i in points]
            for group in groups.values()]
    rhs = [len(group) * demand for (_, demand), group in groups.items()]
    for k in range(1, 4):
        for selection in combinations(range(len(rows)), k):
            mat = [rows[j] for j in selection]
            system = [[Fraction(sum(x*y for x,y in zip(row,col))) for col in mat]
                      + [Fraction(rhs[j])] for row,j in zip(mat,selection)]
            for pivot in range(k):
                diagonal = system[pivot][pivot]
                if not diagonal:
                    break
                basis = [v/diagonal for v in system[pivot]]
                for t in range(k):
                    if t != pivot:
                        scale = system[t][pivot]
                        system[t] = [v-scale*w for v,w in zip(system[t],basis)]
                system[pivot] = basis
            else:
                dual = [r[-1] for r in system]
                if min(dual) < 0:
                    continue
                primal = [sum(x*y for x,y in zip(col,dual)) for col in zip(*mat)]
                if all(sum(p for i,p in zip(points,primal) if outside >> i & 1) >= demand
                       for outside,demand in constraints):
                    return str(sum(p*p for p in primal)), sorted(set(map(str,primal)))
    raise ArithmeticError(mask)


def main():
    counts = walsh_supports()
    gens = [[i ^ 1 for i in range(16)],
            [((i << 1) & 15) | (i >> 3) for i in range(16)],
            [i ^ (2 if i & 1 else 0) for i in range(16)]]
    remaining = set(counts)
    results = []
    while remaining:
        mask = min(remaining)
        copies = orbit(mask, gens, lambda m,gs:(permute(m,g) for g in gs))
        assert copies <= remaining
        remaining.difference_update(copies)
        norm2, weights = certificate(mask)
        results.append(dict(mask=mask, size=mask.bit_count(), orbit=len(copies),
                            sign_pairs=sum(counts[m] for m in copies), b_squared=norm2,
                            optimal_weights=weights))
    assert sum(r['sign_pairs'] for r in results) == 32768
    print(json.dumps(dict(python=platform.python_version(), orbits=results),indent=2))


if __name__ == '__main__':
    main()
