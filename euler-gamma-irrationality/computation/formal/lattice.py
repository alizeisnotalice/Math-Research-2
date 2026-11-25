"""Exact integer kernel, saturation and image-rank certificates."""

import argparse
import json
import sys
from collections import Counter
from functools import cache
from math import comb, gcd, lcm
from pathlib import Path

from flint import fmpz_mat, fmpq_mat

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from intervals import harmonic, certify
from sondow import factors


@cache
def column(n):
    d = lcm(*range(1, 2 * n + 1))
    coeff = [comb(n, i) ** 2 for i in range(n + 1)]
    a = sum(c * int(d * harmonic(n + i)) for i, c in enumerate(coeff))
    prefix = []
    total = 0
    for i in range((n + 1) // 2):
        total += 2 * coeff[i] * int(d * (harmonic(n - i) - harmonic(i)))
        prefix.append(total)
    logs = Counter()
    for k in range(1, n + 1):
        w = prefix[min(k - 1, n - k)]
        for p in factors(n + k):
            logs[p] += w
    return logs, a, d * comb(2 * n, n)


def matrices(n):
    columns = [column(i) for i in range(n, 2 * n + 1)]
    primes = sorted(set().union(*(c[0] for c in columns)))
    E = fmpz_mat([[c[0][p] for c in columns] for p in primes])
    ab = fmpz_mat([[c[j] for c in columns] for j in (1, 2)])
    return E, ab


def saturation(E):
    raw, k = E.nullspace()
    if not k:
        return fmpz_mat(0, E.ncols())
    cols = []
    for j in range(k):
        v = [int(raw[i, j]) for i in range(E.ncols())]
        g = gcd(*v)
        cols.append([x // g for x in v])
    C = fmpz_mat(cols)
    # Columns of C generate an integer lattice with basis H^T.
    h = C.transpose().hnf()
    H = fmpz_mat(h.tolist()[:k]).transpose()
    basis = fmpq_mat(H).inv() * fmpq_mat(C)
    if any(z.q != 1 for row in basis.tolist() for z in row):
        raise ArithmeticError("non-integral saturation")
    K = fmpz_mat([[int(z.p) for z in row] for row in basis.tolist()])
    return K.lll()


def block(n, full=True):
    E, ab = matrices(n)
    K = saturation(E)
    k = K.nrows()
    assert (E * K.transpose()).is_zero()
    assert E.rank() + k == E.ncols()
    if k:
        assert K.transpose().hnf().tolist()[:k] == [
            [int(i == j) for j in range(k)] for i in range(k)
        ]
    image = ab * K.transpose()
    rank = image.rank() if k else 0
    assert rank == fmpz_mat(E.tolist() + ab.tolist()).rank() - E.rank()
    vectors = []
    for j, row in enumerate(K.tolist()):
        v = list(map(int, row))
        a, b = int(image[0, j]), int(image[1, j])
        record = dict(length=sum(map(abs, v)), A=a, B=b)
        if full:
            record.update(v=v, form=certify(b, -a, {}, start=128, limit=512))
        vectors.append(record)
    result = dict(N=n, rank=E.rank(), nullity=k, image_rank=rank, vectors=vectors)
    if k == 2:
        result["image_det"] = int(image.det())
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("N", nargs="+", type=int)
    args = parser.parse_args()
    print(json.dumps([block(n) for n in args.N], indent=2))
