"""Targeted MILP support search beyond diagonal domination and monomial cycles.

The model is necessary-only.  Every retained support has an exact positive
kernel certificate for the two-diagonal monomial-row rule; it is not a
coefficient-level factorization or a counterexample.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import lil_matrix


def load_domination(path):
    spec = importlib.util.spec_from_file_location("domination", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    module.self_check()
    return module


def terms_and_units(I, J):
    edges = defaultdict(list)
    for ii in I:
        for jj in J:
            edges[ii + jj].append((ii, jj))
    edges = dict(sorted(edges.items()))
    units = {("a", 0), ("a", I[-1]), ("b", 0), ("b", J[-1])}
    for pairs in edges.values():
        if len(pairs) == 1:
            ii, jj = pairs[0]
            units.update({("a", ii), ("b", jj)})
    def reduced(edge):
        ii, jj = edge
        return tuple(x for x in (("a", ii), ("b", jj)) if x not in units)
    direct = None
    for k, pairs in edges.items():
        reduced_terms = [reduced(edge) for edge in pairs]
        if len(pairs) > 1 and any(not term for term in reduced_terms):
            direct = {"diagonal": k, "reason": "a forced-unit monomial plus a strictly positive other term"}
            break
    variables = tuple(sorted(({("a", ii) for ii in I} | {("b", jj) for jj in J}) - units))
    return edges, units, variables, reduced, direct


def monomial_rows(I, J):
    edges, units, variables, reduced, direct = terms_and_units(I, J)
    rows = []
    for k, ell in itertools.combinations(edges, 2):
        left = [reduced(edge) for edge in edges[k]]
        right = [reduced(edge) for edge in edges[ell]]
        for mono in sorted(set(left) & set(right)):
            left.remove(mono); right.remove(mono)
        if len(left) == len(right) == 1:
            vector = tuple(left[0].count(var) - right[0].count(var) for var in variables)
            if any(vector):
                rows.append({"diagonals": [k, ell], "left": left[0], "right": right[0], "row": vector})
    return edges, units, variables, direct, rows


def positive_kernel(rows, dimension):
    """Give exact w>0 with Lw=0, proving rowspace(L)∩R_+^d={0}."""
    if dimension == 0:
        return {"status": "no_strict_variables"}
    if not rows:
        return {"status": "exact", "w": [1] * dimension, "reason": "zero row space"}
    L = sp.Matrix([row["row"] for row in rows])
    basis = sp.Matrix.hstack(*L.nullspace())
    if basis.cols == 0:
        return {"status": "none"}
    B = np.array(basis.tolist(), dtype=float)
    # q and delta: maximize delta, Bq>=delta*1, sum(Bq)=1.
    r = basis.cols
    c = np.r_[np.zeros(r), -1.0]
    A_ub = np.c_[-B, np.ones((dimension, 1))]
    A_eq = np.c_[np.sum(B, axis=0, keepdims=True), np.zeros((1, 1))]
    lp = linprog(c, A_ub=A_ub, b_ub=np.zeros(dimension), A_eq=A_eq, b_eq=[1.0],
                 bounds=[(None, None)] * r + [(0, None)], method="highs")
    if not lp.success or lp.x[-1] <= 1e-9:
        return {"status": "none"}
    # A rational q gives an exact kernel vector; use the numerical interior point
    # only to choose it, then verify strict positivity with SymPy rationals.
    q = [sp.Rational(Fraction(float(v)).limit_denominator(10**7)) for v in lp.x[:r]]
    w = basis * sp.Matrix(q)
    if not all(value > 0 for value in w):
        return {"status": "numeric-positive-but-no-rational-reconstruction", "delta": float(lp.x[-1])}
    den = sp.ilcm(*[value.q for value in w])
    wint = [int(den * value) for value in w]
    assert all(value > 0 for value in wint) and L * sp.Matrix(wint) == sp.zeros(L.rows, 1)
    return {"status": "exact", "w": wint, "delta_numeric": float(lp.x[-1])}


def assess(I, J, domination):
    # Standard singleton-unit direct test first.
    edges, units, variables, direct, rows = monomial_rows(I, J)
    if direct:
        return {"status": "direct_unit_budget", "direct": direct}
    cert = domination.find_certificate(I, J)
    if cert is not None:
        domination.verify_certificate(I, J, cert)
        return {"status": "diagonal_domination", "certificate": cert}
    kernel = positive_kernel(rows, len(variables))
    if kernel["status"] == "exact":
        return {"status": "escapes_both_rules", "units": sorted("%s%d" % v for v in units),
                "strict_variables": ["%s%d" % v for v in variables], "monomial_rows": rows,
                "positive_kernel": kernel}
    return {"status": "monomial_rule_not_certified_escape", "kernel": kernel,
            "row_count": len(rows), "strict_variable_count": len(variables)}


def build_and_search(m, n, domination, min_support, limit, seconds, fractional_internal):
    na, nb = m + 1, n + 1
    w0, nw = na + nb, na * nb
    inner_diagonals = [k for k in range(1, m+n) if k not in (m, n)] if fractional_internal else []
    y0 = w0 + nw
    nv = y0 + len(inner_diagonals)
    wij = lambda ii, jj: w0 + ii * nb + jj
    diagonals = {k: [(ii, jj) for ii in range(na) for jj in range(nb) if ii + jj == k]
                 for k in range(m + n + 1)}
    rows, lo, hi = [], [], []
    def add(entries, lower=-np.inf, upper=np.inf):
        rows.append(entries); lo.append(lower); hi.append(upper)
    for ii in range(na):
        for jj in range(nb):
            z = wij(ii, jj)
            add({z: 1, ii: -1}, upper=0)
            add({z: 1, na + jj: -1}, upper=0)
            add({z: 1, ii: -1, na + jj: -1}, lower=-1)
    # Only the endpoint-adjacent diagonals must be singleton; other singleton
    # diagonals are deliberately allowed and handled by exact unit closure later.
    for k in (m, n):
        add({wij(ii, jj): 1 for ii, jj in diagonals[k]}, lower=1, upper=1)
    if fractional_internal:
        for offset, k in enumerate(inner_diagonals):
            yy = y0 + offset
            terms = {wij(ii, jj): 1 for ii, jj in diagonals[k]}
            add({**terms, yy: -len(diagonals[k])}, upper=0)
            add({**terms, yy: -2}, lower=0)
    add({wij(t, t): 1 for t in range(1, m + 1)}, lower=1)  # positive intersection
    add({q: 1 for q in range(na + nb)}, lower=min_support)
    lower = np.zeros(nv); upper = np.ones(nv)
    for q in (0, m, na, na + n): lower[q] = upper[q] = 1
    objective = np.zeros(nv); objective[:na + nb] = 1
    deadline = time.monotonic() + seconds
    records, counts = [], defaultdict(int)
    for round_no in range(limit):
        M = lil_matrix((len(rows), nv), dtype=float)
        for rr, terms in enumerate(rows):
            for col, value in terms.items(): M[rr, col] = value
        result = milp(c=objective, integrality=np.ones(nv), bounds=Bounds(lower, upper),
                      constraints=LinearConstraint(M.tocsr(), lo, hi),
                      options={"time_limit": max(1.0, deadline-time.monotonic()), "mip_rel_gap": 0})
        if result.x is None:
            records.append({"round": round_no, "solver_status": int(result.status), "solver_message": result.message})
            break
        bits = np.rint(result.x[:na+nb]).astype(int)
        I = tuple(q for q in range(na) if bits[q])
        J = tuple(q for q in range(nb) if bits[na+q])
        verdict = assess(I, J, domination)
        counts[verdict["status"]] += 1
        records.append({"round": round_no, "I": I, "J": J, "support_size": len(I)+len(J), "verdict": verdict})
        selected = [q for q, bit in enumerate(bits) if bit]
        unselected = [q for q, bit in enumerate(bits) if not bit]
        add({**{q: 1 for q in selected}, **{q: -1 for q in unselected}}, upper=len(selected)-1)
        if time.monotonic() >= deadline: break
    return {"m": m, "n": n, "min_support": min_support, "limit": limit,
            "model": ("endpoint singleton plus internal r=0-or-at-least-2 supports" if fractional_internal
                      else "endpoint singleton plus unrestricted internal singleton supports"),
            "counts": dict(counts), "records": records, "elapsed_seconds": seconds-max(0, deadline-time.monotonic())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, required=True); ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--min-support", type=int, default=11); ap.add_argument("--limit", type=int, default=80)
    ap.add_argument("--seconds", type=float, default=1800); ap.add_argument("--domination", type=Path, required=True)
    ap.add_argument("--fractional-internal", action="store_true")
    args = ap.parse_args()
    if not (2 < args.m < args.n and args.limit > 0 and args.seconds > 0): ap.error("invalid range")
    domination = load_domination(args.domination)
    data = {"args": vars(args), "python": sys.version, "scipy": __import__("scipy").__version__,
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "result": build_and_search(args.m, args.n, domination, args.min_support, args.limit, args.seconds,
                                         args.fractional_internal)}
    print(json.dumps(data, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__": main()
