#!/usr/bin/env python3
"""Exact and small CPU audits for the tau-stratified original-law particle model.

The reference law is a uniformly random descending vertex order.  This script
does not estimate the original law from a fixed-prefix completion kernel.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from math import factorial
from pathlib import Path
import json
import random
import sys


def graph(n: int):
    N = n * n
    adj = [0] * N
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    adj[v] |= 1 << (yy * n + xx)
    return adj


def boundary(mask: int, adj: list[int]) -> list[int]:
    return [v for v in range(len(adj)) if not (mask >> v & 1) and adj[v] & mask]


def birth_set(mask: int, adj: list[int]) -> list[int]:
    return [v for v in range(len(adj)) if not (mask >> v & 1) and not (adj[v] & mask)]


def distance(n: int, u: int, v: int) -> int:
    return abs(u % n - v % n) + abs(u // n - v // n)


def squared_euclidean(n: int, u: int, v: int) -> int:
    """The Q_n normalization used in the tail bound, not Manhattan squared."""
    dx = u % n - v % n
    dy = u // n - v // n
    return dx * dx + dy * dy


def exact_permutations(n: int):
    """Return integer counts indexed by tau (one-based descending rank)."""
    adj, N = graph(n), n * n
    count = [0] * (N + 1)
    dist_sum = [0] * (N + 1)
    for order in permutations(range(N)):
        mask = 0
        births = 0
        roots: list[int] = []
        tau = 0
        for t, v in enumerate(order, 1):
            if not (adj[v] & mask):
                births += 1
                roots.append(v)
                if births == 2:
                    tau = t
            mask |= 1 << v
        if births == 2:
            count[tau] += 1
            dist_sum[tau] += squared_euclidean(n, roots[0], roots[1])
    return count, dist_sum


def proposal_tree(n: int, tau: int):
    """Enumerate q_tau paths and sum q_tau*W_tau exactly.

    Every final constrained order has q*W=1/N!, so this independently audits
    the Feynman--Kac potential convention and the tau event.
    """
    adj, N = graph(n), n * n
    z = Fraction(0)
    m = Fraction(0)

    # Carry both births in the recursion explicitly.  The constrained path
    # tree is intentionally distinct from the all-permutation enumeration.
    def rec(t: int, mask: int, root: int, second_root: int, w: Fraction):
        nonlocal z, m
        if t == N:
            z += w
            m += w * squared_euclidean(n, root, second_root)
            return
        rank = t + 1
        if rank == 1:
            allowed = list(range(N))
        elif rank == tau:
            allowed = birth_set(mask, adj)
        else:
            allowed = boundary(mask, adj)
        if not allowed:
            return
        g = Fraction(len(allowed), N - t)
        for v in allowed:
            rec(rank, mask | (1 << v), v if rank == 1 else root,
                v if rank == tau else second_root, w * g * Fraction(1, len(allowed)))

    rec(0, 0, -1, -1, Fraction(1))
    return z, m


MASK64 = (1 << 64) - 1


def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    return (x ^ (x >> 31)) & MASK64


def smc_once(n: int, tau: int, particles: int, seed: int):
    """Bootstrap particle filter for a fixed tau stratum.

    The particle potential before rank r is |A_r|/(N-r+1), where A_r is the
    forced allowed set.  Multinomial resampling is used because its standard
    first-moment identity yields unbiased normalizing and terminal unnormalized
    estimators.  Child RNG keys are regenerated from (seed, rank, child), not
    copied from a parent.
    """
    adj, N = graph(n), n * n
    states = [(0, -1, -1)] * particles  # revealed mask, first root, second root
    zhat = 1.0
    min_unique_ancestors = particles
    for t in range(N):
        rank = t + 1
        allowed: list[list[int]] = []
        g: list[float] = []
        for mask, _, _ in states:
            if rank == 1:
                a = list(range(N))
            elif rank == tau:
                a = birth_set(mask, adj)
            else:
                a = boundary(mask, adj)
            allowed.append(a)
            g.append(len(a) / (N - t))
        mean_g = sum(g) / particles
        if mean_g == 0.0:
            return 0.0, 0.0, 0
        # At the final rank retain the weighted particle average instead of
        # resampling once more.  This is the direct terminal Feynman--Kac
        # estimator z_{N-1} P^{-1} sum_i G_N(i) f(X_i), and removes one
        # unnecessary resampling noise layer.
        if t == N - 1:
            weighted_f = 0.0
            for child, ((mask, root1, root2), a, potential) in enumerate(zip(states, allowed, g)):
                if not a:
                    continue
                rng = random.Random(splitmix64(seed ^ (rank << 40) ^ (child * 0x9E37)))
                v = a[rng.randrange(len(a))]
                r1 = v if rank == 1 else root1
                r2 = v if rank == tau else root2
                weighted_f += potential * squared_euclidean(n, r1, r2)
            return zhat * mean_g, zhat * weighted_f / particles, min_unique_ancestors
        zhat *= mean_g
        cdf = []
        acc = 0.0
        for x in g:
            acc += x / (particles * mean_g)
            cdf.append(acc)
        parent = []
        for child in range(particles):
            rng = random.Random(splitmix64(seed ^ (rank << 32) ^ child))
            u = rng.random()
            lo, hi = 0, particles - 1
            while lo < hi:
                mid = (lo + hi) // 2
                if u <= cdf[mid]:
                    hi = mid
                else:
                    lo = mid + 1
            parent.append(lo)
        min_unique_ancestors = min(min_unique_ancestors, len(set(parent)))
        next_states = []
        for child, p in enumerate(parent):
            mask, root1, root2 = states[p]
            # This key is intentionally not inherited from the selected parent.
            rng = random.Random(splitmix64(seed ^ (rank << 40) ^ (child * 0x9E37)))
            v = allowed[p][rng.randrange(len(allowed[p]))]
            next_states.append((mask | (1 << v), v if rank == 1 else root1,
                                v if rank == tau else root2))
        states = next_states
    mean_f = sum(squared_euclidean(n, r1, r2) for _, r1, r2 in states) / particles
    return zhat, zhat * mean_f, min_unique_ancestors


def main(out_path: Path):
    result = {"implementation": "Python 3 standard library", "exact": {}, "tree_audit": {}, "smc": {}}
    for n in (2, 3):
        count, dsum = exact_permutations(n)
        N = n * n
        exact_z = [Fraction(x, factorial(N)) for x in count]
        exact_m = [Fraction(x, factorial(N)) for x in dsum]
        assert sum(count) == (8 if n == 2 else 143808)
        result["exact"][str(n)] = {
            "A2": sum(count),
            "A2_by_tau": {str(t): count[t] for t in range(2, N + 1) if count[t]},
            "dist2_sum_by_tau": {str(t): dsum[t] for t in range(2, N + 1) if dsum[t]},
            "Q": str(sum(exact_m) / (sum(exact_z) * n * n)),
        }
        for tau in range(2, N + 1):
            if count[tau]:
                z, m = proposal_tree(n, tau)
                assert z == exact_z[tau]
                assert m == exact_m[tau]
        result["tree_audit"][str(n)] = "all nonzero tau strata exactly agree"
        reps, P = 128, 256
        vals = {tau: [] for tau in range(2, N + 1) if count[tau]}
        for rep in range(reps):
            for tau in vals:
                vals[tau].append(smc_once(n, tau, P, 0xA0220000 + 10000 * n + 100 * tau + rep))
        zbar = sum(sum(x[0] for x in rows) / reps for rows in vals.values())
        mbar = sum(sum(x[1] for x in rows) / reps for rows in vals.values())
        result["smc"][str(n)] = {
            "particles": P,
            "replicates": reps,
            "Z_mean": zbar,
            "M_mean": mbar,
            "Q_ratio": mbar / (zbar * n * n),
            "min_unique_ancestors": min(x[2] for rows in vals.values() for x in rows),
            "note": "finite repeated-ratio diagnostic only; Z and M are the unnormalized estimands",
        }
    out_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("tau_smc_audit.json")
    main(target)
