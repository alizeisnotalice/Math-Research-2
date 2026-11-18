"""Exact integer audits for block-rank and uniform-affine cut claims."""
from fractions import Fraction
import itertools
import json
import random
import sys

from exact_rank6 import points, rank, subspaces

SEED = 8527


def form(size, rng):
    rows = [0] * size
    for i, j in itertools.combinations(range(size), 2):
        if rng.getrandbits(1):
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return tuple(rows)


def schur_checks(rng):
    checked = 0
    for a, n, j in itertools.product(range(6), repeat=3):
        for _ in range(12):
            size = a + n + j
            M = list(form(size, rng))
            whole_q = (1 << (a + n)) - 1
            only_a = (1 << a) - 1
            only_j = ((1 << size) - 1) ^ whole_q
            for i in range(a, a + n):
                M[i] &= only_j
            for i in range(a):
                M[i] &= ~(whole_q ^ only_a)
            r = rank(M)
            rho = rank(row & only_a for row in M[:a])
            k = rank(row >> (a + n) for row in M[a:a+n])
            p = rank(row & whole_q for row in M)
            assert p >= rho + k and r >= rho + 2*k
            reserve = (1 << (r//2)) - 1
            inside = (1 << (r//2)) * (1 - Fraction(1, 1 << p))
            assert reserve >= ((1 << (rho//2))-1)+((1<<k)-1)
            phi = (1 << (rho//2)) - Fraction(1, 1 << (rho//2))
            assert inside >= phi + ((1<<k)-1)
            checked += 1
    return checked


def expand_pencil(generators):
    result = [generators[0]]
    for delta in generators[1:]:
        result.extend(tuple(x ^ y for x, y in zip(old,delta)) for old in result[:])
    return result


def uniform_cuts(rng):
    examples = 0
    for v in (2, 4, 6):
        for h in (4, 6, 8, 12, 16):
            for _ in range(12):
                pencil = expand_pencil([form(h, rng) for _ in range(v+1)])
                if rank(row for matrix in pencil for row in matrix) < h:
                    continue
                r = [rank(matrix) for matrix in pencil]
                for t in range(h+1):
                    mask = (1<<t)-1
                    p = [rank(row & mask for row in matrix) for matrix in pencil]
                    # Coordinate cuts on the old frequency; include both U=V and U proper.
                    for u in range(v+1):
                        if u+t == 0:
                            continue
                        g = Fraction(0)
                        for e, (polar, projected) in enumerate(zip(r,p)):
                            inside = not (e & ((1<<u)-1))
                            weight = 1 << (polar//2)
                            g += weight*(1-Fraction(1,1<<projected)) if inside else weight
                        assert (v+2)*g >= (1<<v)*(u+t+1)
                        examples += 1
    return examples


if __name__ == '__main__':
    rng = random.Random(SEED)
    print(json.dumps({"python":sys.version.split()[0], "seed":SEED,
                      "integer_schur_checks":schur_checks(rng),
                      "uniform_affine_projected_cuts":uniform_cuts(rng)}, indent=2))
