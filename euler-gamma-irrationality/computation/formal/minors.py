"""Modular full-row-rank witnesses for E augmented by the two arithmetic rows."""

import json
from math import isqrt
from flint import fmpz_mat, nmod_mat
from lattice import matrices


def prime(x):
    return x >= 2 and all(x % i for i in range(2, isqrt(x) + 1))


def witness(N):
    E, AB = matrices(N)
    rows = E.tolist() + AB.tolist()
    q = 4 * N + 1
    while not prime(q):
        q += 1
    mod = nmod_mat([[int(x) % q for x in r] for r in rows], q)
    echelon, rank = mod.rref()
    pivots = [next(i for i, x in enumerate(r) if x) for r in echelon.tolist()[:rank]]
    if rank != len(rows):
        return dict(N=N, q=q, rank=rank, rows=len(rows), pivots=pivots)
    sub = nmod_mat([[int(r[i]) % q for i in pivots] for r in rows], q)
    det = int(sub.det())
    assert det != 0
    return dict(N=N, q=q, rank=rank, rows=len(rows), pivots=pivots, det_mod_q=det)


if __name__ == "__main__":
    print(
        json.dumps(
            [witness(N) for N in (32, 33, 35, 38, 40, 48, 64, 96, 128)], indent=2
        )
    )
