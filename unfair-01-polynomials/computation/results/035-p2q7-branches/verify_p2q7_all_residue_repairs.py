#!/usr/bin/env python3
"""Independent coverage and status verifier for all-residue p2q7 audit."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

P, Q = 2, 7
TS = {(0, 1, 4), (0, 1), (0, 3, 4), (0, 3, 6), (0, 3), (0, 4), (0, 6)}


def G_of(T): return tuple(sorted(set(T) | {x + Q for x in T}))


def opts_of(G):
    ans = []
    for residue in range(P):
        members = [x for x in G if x % P == residue]
        if members:
            j0 = min(members)
            ans.append([(x, Q * (x - j0) // P) for x in members if x > j0])
    return ans


def status(T, selected, ordinary_j):
    G = G_of(T)
    D = {P, ordinary_j, *(offset for _, offset in selected)}
    terms = {}
    for d in sorted(D):
        for g in G:
            if g < d:
                terms.setdefault(d - g, []).append((d, g))
    overlaps = any(len(v) > 1 for v in terms.values())
    offsets = sorted(terms)
    pp, qp = set(), set()
    for i, b in enumerate(offsets):
        for a in offsets[:i]:
            delta = b-a
            if delta % P == 0: pp.add(delta // P)
            if delta % Q == 0: qp.add(delta // Q)
    return overlaps or bool(pp & qp), terms


def main() -> None:
    source = Path("p2q7-all-residue-search.json")
    data = json.loads(source.read_text())
    assert data["input"]["p"]==P and data["input"]["q"]==Q
    assert {tuple(x) for x in data["input"]["T_patterns"]} == TS
    expected = set()
    for T in TS:
        G = G_of(T)
        options = opts_of(G)
        assert all(options)
        for selected in itertools.product(*options):
            for j in G:
                if j:
                    expected.add((T, selected, j))
    seen = set(); closed = opened = 0
    for row in data["branches"]:
        T = tuple(row["T"])
        selected = tuple((x["G_index"], x["F_offset"]) for x in row["residue_repair_choices"])
        key = (T, selected, row["ordinary_repair_j"])
        assert key not in seen; seen.add(key)
        is_closed, terms = status(*key)
        reported = {x["offset"]: [tuple(y) for y in x["positive_terms"]] for x in row["below_r_positive_terms_by_H_offset"]}
        assert reported == terms
        assert (row["status"] == "closed_by_forced_positive_collision") == is_closed
        closed += is_closed; opened += not is_closed
    assert seen == expected
    assert data["summary"]["branches"] == len(seen)
    assert data["summary"]["closed"] == closed and data["summary"]["open"] == opened
    receipt = {"status": "passed", "branches": len(seen), "closed": closed, "open": opened,
               "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
               "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path("p2q7-all-residue-verifier-receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__": main()
