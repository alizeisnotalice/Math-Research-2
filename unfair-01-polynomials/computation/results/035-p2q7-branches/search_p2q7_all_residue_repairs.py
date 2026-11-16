#!/usr/bin/env python3
"""Complete p=2,q=7 residue-minimum repair branching for the specified G supports.

For each occupied residue class mod p, A at its least G index has a negative
f_r term and therefore requires at least one positive F coefficient at
offset q*(j-j0)/p, for a later G index j in that class.  The search takes one
such necessary repair per class and one ordinary h_r repair.  It uses only
positive products below r; every branch is conditional on r being at least
its selected F offsets, as required for that F coefficient to exist.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path


P, Q = 2, 7
TS = ((0, 1, 4), (0, 1), (0, 3, 4), (0, 3, 6), (0, 3), (0, 4), (0, 6))


def G_of(T: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(set(T) | {x + Q for x in T}))


def residue_options(G: tuple[int, ...]) -> list[dict]:
    out = []
    for residue in range(P):
        members = [j for j in G if j % P == residue]
        if not members:
            continue
        j0 = min(members)
        later = [j for j in members if j > j0]
        out.append({
            "residue": residue,
            "minimum_G_index": j0,
            "choices": [{"G_index": j, "F_offset": Q * (j - j0) // P} for j in later],
        })
    return out


def c_collisions(offsets: list[int]) -> list[dict]:
    pp: dict[int, list[tuple[int, int]]] = {}
    qp: dict[int, list[tuple[int, int]]] = {}
    for i, b in enumerate(offsets):
        for a in offsets[:i]:
            delta = b - a
            if delta % P == 0:
                pp.setdefault(delta // P, []).append((a, b))
            if delta % Q == 0:
                qp.setdefault(delta // Q, []).append((a, b))
    result = []
    for k in sorted(set(pp) & set(qp)):
        for qpair in qp[k]:
            for ppair in pp[k]:
                a, b = qpair
                c, d = ppair
                result.append({
                    "k": k, "q_offset_pair": [a, b], "p_offset_pair": [c, d],
                    "identity": f"2*(r-{a})+7*(r-{d}) = 2*(r-{b})+7*(r-{c})",
                })
    return result


def analyze(T: tuple[int, ...], selected: tuple[dict, ...], ordinary_j: int) -> dict:
    G = G_of(T)
    D = tuple(sorted({P, ordinary_j, *(x["F_offset"] for x in selected)}))
    terms: dict[int, list[tuple[int, int]]] = {}
    for d in D:
        for g in G:
            if g < d:
                terms.setdefault(d - g, []).append((d, g))
    overlaps = [{"offset": o, "positive_terms": v} for o, v in sorted(terms.items()) if len(v) >= 2]
    cs = c_collisions(sorted(terms))
    witnesses = []
    if overlaps:
        witnesses.append({"kind": "binary_H_overlap", "details": overlaps})
    if cs:
        witnesses.append({"kind": "binary_C_collision", "details": cs})
    return {
        "T": list(T), "G_support": list(G), "residue_repair_choices": list(selected),
        "ordinary_repair_j": ordinary_j, "positive_F_offsets_D": list(D),
        "below_r_positive_terms_by_H_offset": [
            {"offset": o, "positive_terms": v} for o, v in sorted(terms.items())
        ],
        "status": "closed_by_forced_positive_collision" if witnesses else "open_under_this_rule_set",
        "witnesses": witnesses,
    }


def main() -> None:
    branches = []
    immediate_singleton = []
    old42_coverage = 0
    for T in TS:
        G = G_of(T)
        opts = residue_options(G)
        if any(not x["choices"] for x in opts):
            immediate_singleton.append(list(T))
            continue
        for selected in itertools.product(*(x["choices"] for x in opts)):
            for j in G:
                if j:
                    branches.append(analyze(T, selected, j))
        # The former 42-branch audit chose only an even A_2r repair and j.
        old42_coverage += sum(1 for v in G if v and v % 2 == 0) * (len(G) - 1)
    summary = {
        "T_patterns": len(TS), "immediate_singleton_T": immediate_singleton,
        "branches": len(branches), "closed": sum(x["status"].startswith("closed") for x in branches),
        "open": sum(x["status"].startswith("open") for x in branches),
        "all_branches_closed": all(x["status"].startswith("closed") for x in branches),
        "previous_even_only_branch_count": old42_coverage,
    }
    result = {
        "schema": "p2q7-all-residue-repairs-v1",
        "input": {"p": P, "q": Q, "T_patterns": [list(x) for x in TS]},
        "scope": "Every branch chooses one necessary later-G repair for each nonempty mod-2 class plus one necessary h_r repair. It is not a coefficient solver.",
        "summary": summary, "branches": branches,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    Path("p2q7-all-residue-search.json").write_text(encoded)
    print(json.dumps(summary, sort_keys=True))
    print("output_sha256=" + hashlib.sha256(encoded.encode()).hexdigest())
    print("source_sha256=" + hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
