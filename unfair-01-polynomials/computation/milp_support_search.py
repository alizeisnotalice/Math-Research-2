"""MILP search for support patterns outside singleton and domination filters.

This searches only a combinatorial necessary-condition model.  A surviving pair
is not a polynomial factorization, nor a counterexample to the conjecture.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix


def load_domination(path: Path):
    spec = importlib.util.spec_from_file_location("diagonal_domination", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    module.self_check()
    return module


def solve_pair(m: int, n: int, domination, limit: int, time_limit: float):
    # Variables: A_i, B_j, W_ij=A_i AND B_j, and binary occupancy Y_k for
    # every non-corner diagonal.  r_k is 0 if Y_k=0, and >=2 if Y_k=1.
    na, nb = m + 1, n + 1
    w0 = na + nb
    nw = na * nb
    y_diag = [k for k in range(1, m + n) if k not in (m, n)]
    y0 = w0 + nw
    nv = y0 + len(y_diag)
    w_index = lambda ii, jj: w0 + ii * nb + jj
    y_index = {k: y0 + q for q, k in enumerate(y_diag)}
    diag_edges = {k: [(ii, jj) for ii in range(na) for jj in range(nb) if ii + jj == k]
                  for k in range(m + n + 1)}

    rows, lo, hi = [], [], []
    def add(entries, lower=-np.inf, upper=np.inf):
        rows.append(entries); lo.append(lower); hi.append(upper)

    # AND linearization.
    for ii in range(na):
        for jj in range(nb):
            wij = w_index(ii, jj)
            add({wij: 1, ii: -1}, upper=0)
            add({wij: 1, na + jj: -1}, upper=0)
            add({wij: 1, ii: -1, na + jj: -1}, lower=-1)
    # Endpoint singleton diagonals.
    for k in (m, n):
        add({w_index(ii, jj): 1 for ii, jj in diag_edges[k]}, lower=1, upper=1)
    # Every other inner diagonal is empty or has at least two support pairs.
    for k in y_diag:
        terms = {w_index(ii, jj): 1 for ii, jj in diag_edges[k]}
        yy = y_index[k]
        # r_k <= M*y_k and r_k >= 2*y_k.
        add({**terms, yy: -len(diag_edges[k])}, upper=0)
        add({**terms, yy: -2}, lower=0)
    # Require at least one positive common index.  W_tt encodes A_t=B_t=1.
    common = {w_index(t, t): 1 for t in range(1, min(m, n) + 1)}
    add(common, lower=1)

    base_rows = list(rows); base_lo = list(lo); base_hi = list(hi)
    lower = np.zeros(nv); upper = np.ones(nv)
    for idx in (0, m, na, na + n):
        lower[idx] = upper[idx] = 1
    integrality = np.ones(nv)
    objective = np.zeros(nv)
    objective[:na + nb] = 1
    candidates = []
    excluded = 0
    solve_log = []
    deadline = time.monotonic() + time_limit
    for round_no in range(limit):
        remaining = max(1.0, deadline - time.monotonic())
        if remaining <= 1.0 and round_no:
            break
        M = lil_matrix((len(rows), nv), dtype=float)
        for rr, entries in enumerate(rows):
            for col, val in entries.items():
                M[rr, col] = val
        result = milp(c=objective, integrality=integrality, bounds=Bounds(lower, upper),
                      constraints=LinearConstraint(M.tocsr(), np.asarray(lo), np.asarray(hi)),
                      options={"time_limit": remaining, "mip_rel_gap": 0.0})
        solve_log.append({"round": round_no, "status": int(result.status), "message": result.message,
                          "fun": None if result.fun is None else float(result.fun)})
        if result.x is None:
            break
        bits = np.rint(result.x[:na + nb]).astype(int)
        I = tuple(k for k in range(na) if bits[k])
        J = tuple(k for k in range(nb) if bits[na + k])
        cert = domination.find_certificate(I, J)
        if cert is None:
            candidates.append({"I": I, "J": J, "support_size": len(I)+len(J),
                               "domination_status": "unknown"})
            if len(candidates) >= limit:
                break
        else:
            domination.verify_certificate(I, J, cert)
            excluded += 1
        # No-good: differ in at least one A/B support bit.
        selected = [q for q, bit in enumerate(bits) if bit]
        unselected = [q for q, bit in enumerate(bits) if not bit]
        add({**{q: 1 for q in selected}, **{q: -1 for q in unselected}}, upper=len(selected)-1)
    return {"m": m, "n": n, "model": "MILP support necessary conditions", "candidates": candidates,
            "domination_excluded": excluded, "solver_log": solve_log,
            "elapsed_seconds": time_limit - max(0.0, deadline-time.monotonic())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=7)
    ap.add_argument("--n", type=int, default=17)
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--seconds", type=float, default=900)
    ap.add_argument("--domination", type=Path, required=True)
    args = ap.parse_args()
    if not (1 < args.m < args.n and args.limit > 0 and args.seconds > 0):
        ap.error("need 1<m<n, positive limit and seconds")
    domination = load_domination(args.domination)
    output = {"args": vars(args), "python": sys.version,
              "scipy": __import__("scipy").__version__,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "result": solve_pair(args.m, args.n, domination, args.limit, args.seconds)}
    print(json.dumps(output, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
