#!/usr/bin/env python3
"""Exact proposal-tree audit for the original K=1/K=2 peak laws.

The orders here are *descending* label orders.  A birth is a vertex with no
previously revealed neighbour, hence it is a peak in the final labelling.
The proposal is deliberately the one intended for later SIS/GPU use:
  b=0: choose any vertex;
  b=1: choose any unrevealed vertex;
  b=2: choose only the current boundary.
Every K<=2 order has positive proposal probability and no other order does.
All arithmetic below is exact rational arithmetic.
"""
from collections import Counter
from fractions import Fraction
import argparse
import json


def graph(n):
    N = n * n
    adj = [0] * N
    xy = []
    for y in range(n):
        for x in range(n):
            v = y * n + x
            xy.append((x, y))
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                X, Y = x + dx, y + dy
                if 0 <= X < n and 0 <= Y < n:
                    adj[v] |= 1 << (Y * n + X)
    return xy, adj


def vertices(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def audit(n):
    xy, adj = graph(n)
    N = n * n
    full = (1 << N) - 1
    terminal = Counter()
    strata = Counter()
    sq_by_stratum = Counter()
    proposal_mass = Fraction(0)
    weighted_A1 = Fraction(0)
    weighted_A2 = Fraction(0)
    weighted_euclidean_sq = Fraction(0)
    weighted_manhattan_sq = Fraction(0)
    leaves = 0

    def rec(S, peaks, tau, prob, weight):
        nonlocal proposal_mass, weighted_A1, weighted_A2, weighted_euclidean_sq, weighted_manhattan_sq, leaves
        if S == full:
            leaves += 1
            assert prob * weight == 1
            proposal_mass += prob
            k = len(peaks)
            terminal[k] += 1
            if k == 1:
                weighted_A1 += prob * weight
            elif k == 2:
                u, v = peaks
                dx = abs(xy[u][0] - xy[v][0])
                dy = abs(xy[u][1] - xy[v][1])
                rho = dx + dy
                e2 = dx * dx + dy * dy
                weighted_A2 += prob * weight
                weighted_euclidean_sq += prob * weight * e2
                weighted_manhattan_sq += prob * weight * rho * rho
                strata[tau] += 1
                sq_by_stratum[tau] += e2
            else:
                raise AssertionError(k)
            return
        if not peaks:
            allowed = full ^ S
        elif len(peaks) == 1:
            allowed = full ^ S
        else:
            # A non-boundary choice would be an impermissible third peak.
            allowed = 0
            for u in vertices(S):
                allowed |= adj[u]
            allowed &= full ^ S
        a = allowed.bit_count()
        assert a > 0
        for v in vertices(allowed):
            born = not (adj[v] & S)
            new_peaks = peaks + ((v,) if born else ())
            assert len(new_peaks) <= 2
            new_tau = tau
            if len(peaks) == 1 and born:
                new_tau = S.bit_count() + 1  # descending rank of peak two
            rec(S | (1 << v), new_peaks, new_tau, prob / a, weight * a)

    rec(0, (), None, Fraction(1), 1)
    assert proposal_mass == 1
    assert weighted_A1.denominator == weighted_A2.denominator == weighted_euclidean_sq.denominator == weighted_manhattan_sq.denominator == 1
    result = {
        "n": n,
        "proposal": "b0:any, b1:any, b2:boundary",
        "proposal_terminal_mass": str(proposal_mass),
        "proposal_tree_leaves": leaves,
        "unweighted_terminal_orders": {str(k): terminal[k] for k in sorted(terminal)},
        "A1_weighted": int(weighted_A1),
        "A2_weighted": int(weighted_A2),
        "euclidean_squared_total_weighted": int(weighted_euclidean_sq),
        "manhattan_squared_total_weighted": int(weighted_manhattan_sq),
        "Q_euclidean": str(Fraction(weighted_euclidean_sq, N * weighted_A2)),
        "Q_manhattan": str(Fraction(weighted_manhattan_sq, N * weighted_A2)),
        "second_birth_rank_counts": dict(sorted(strata.items())),
        "second_birth_rank_squared_totals": dict(sorted(sq_by_stratum.items())),
    }
    return result


def audit_forced_stratum(n, tau):
    """Exact tree for the proposal conditioned to put birth two at rank tau."""
    xy, adj = graph(n)
    N = n * n
    full = (1 << N) - 1
    probability_mass = Fraction(0)
    A = Fraction(0)
    D = Fraction(0)
    leaves = 0

    def boundary(S):
        ans = 0
        for u in vertices(S):
            ans |= adj[u]
        return ans & (full ^ S)

    def rec(S, peaks, prob, weight):
        nonlocal probability_mass, A, D, leaves
        rank = S.bit_count() + 1
        if S == full:
            assert len(peaks) == 2 and prob * weight == 1
            probability_mass += prob
            leaves += 1
            u, v = peaks
            e2 = (xy[u][0] - xy[v][0]) ** 2 + (xy[u][1] - xy[v][1]) ** 2
            A += prob * weight
            D += prob * weight * e2
            return
        if rank == 1:
            allowed = full ^ S
        elif rank < tau:
            allowed = boundary(S)
        elif rank == tau:
            allowed = (full ^ S) ^ boundary(S)
        else:
            allowed = boundary(S)
        a = allowed.bit_count()
        if not a:
            return
        for v in vertices(allowed):
            born = not (adj[v] & S)
            new_peaks = peaks + ((v,) if born else ())
            # By construction every legal branch has exactly one/two births here.
            assert len(new_peaks) <= 2
            rec(S | (1 << v), new_peaks, prob / a, weight * a)

    rec(0, (), Fraction(1), 1)
    # A forced prefix can dead-end before tau; it is a zero-contribution
    # outcome, so the successful-leaf mass may be strictly between 0 and 1.
    assert 0 <= probability_mass <= 1
    assert A.denominator == D.denominator == 1
    return {"proposal_mass": str(probability_mass), "leaves": leaves,
            "A2": int(A), "euclidean_squared_total": int(D)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int, choices=(2, 3))
    args = parser.parse_args()
    out = audit(args.n)
    expected = {
        2: (16, 8, Fraction(1, 2)),
        3: (29568, 143808, Fraction(1100, 2247)),
    }
    A1, A2, Q = expected[args.n]
    assert (out["A1_weighted"], out["A2_weighted"], Fraction(out["Q_euclidean"])) == (A1, A2, Q)
    forced = {tau: audit_forced_stratum(args.n, tau) for tau in range(2, args.n * args.n + 1)}
    assert {tau: row["A2"] for tau, row in forced.items() if row["A2"]} == out["second_birth_rank_counts"]
    assert {tau: row["euclidean_squared_total"] for tau, row in forced.items() if row["A2"]} == out["second_birth_rank_squared_totals"]
    out["forced_second_birth_strata"] = forced
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
