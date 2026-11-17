#!/usr/bin/env python3
"""Independent subset-DP verification of the finite seven-vertex R45 lemma.

For each ordered non-edge (r,p), the state is an already exposed vertex mask
and the first full-rank neighbour of p.  An extension is legal exactly when
the new vertex is adjacent to the mask.  Thus every DP path is one admissible
ordered prefix, with no permutation enumeration used in the full 5184-graph
scan.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import json
import platform
import sys
import time


HERE = Path(__file__).resolve().parent
N = 7


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def adjacency(n: int, edges: tuple[tuple[int, int], ...]) -> list[int]:
    result = [0] * n
    for x, y in edges:
        assert 0 <= x < n and 0 <= y < n and x != y
        result[x] |= 1 << y
        result[y] |= 1 << x
    return result


def h_by_subset_dp(adj: list[int]) -> tuple[list[int], dict[str, int]]:
    """Count h_L by the specified connected-prefix subset DP using integers."""
    n, allmask = len(adj), (1 << len(adj)) - 1
    by_size = [[mask for mask in range(1 << n) if mask.bit_count() == k] for k in range(n + 1)]
    neighbor_union = [0] * (1 << n)
    for mask in range(1, 1 << n):
        bit = mask & -mask
        neighbor_union[mask] = neighbor_union[mask ^ bit] | adj[bit.bit_length() - 1]
    h, starts, states, transitions = [0] * n, 0, 0, 0
    for r in range(n):
        for p in range(n):
            if r == p or (adj[r] >> p) & 1:
                continue
            starts += 1
            initial = (1 << r) | (1 << p)
            # slots 0,...,n-2 encode L=0 (not yet seen) and L=1,...,n-2.
            dp: list[list[int] | None] = [None] * (1 << n)
            dp[initial] = [1] + [0] * (n - 2)
            for rank in range(2, n):
                for mask in by_size[rank]:
                    row = dp[mask]
                    if row is None:
                        continue
                    states += sum(x != 0 for x in row)
                    choices = neighbor_union[mask] & ~mask & allmask
                    while choices:
                        bit = choices & -choices
                        choices ^= bit
                        newmask = mask | bit
                        target = dp[newmask]
                        if target is None:
                            target = [0] * (n - 1)
                            dp[newmask] = target
                        p_is_adjacent = bool(adj[p] & bit)
                        for ell, count in enumerate(row):
                            if not count:
                                continue
                            newell = ell if ell else ((rank + 1) - 2 if p_is_adjacent else 0)
                            target[newell] += count
                            transitions += 1
            final = dp[allmask]
            if final is not None:
                for ell in range(1, n - 1):
                    h[ell] += final[ell]
    return h, {"ordered_nonedge_starts": starts, "nonzero_dp_cells": states, "integer_transitions": transitions}


def h_by_direct_words(adj: list[int]) -> list[int]:
    """Small reference only: direct words for two cross-check graphs, not the scan."""
    n, result = len(adj), [0] * len(adj)
    for word in permutations(range(n)):
        r, p = word[:2]
        if (adj[r] >> p) & 1:
            continue
        mask, ell, valid = (1 << r) | (1 << p), 0, True
        for zero_index, v in enumerate(word[2:], start=2):
            bit = 1 << v
            if not adj[v] & mask:
                valid = False
                break
            if ell == 0 and adj[p] & bit:
                ell = (zero_index + 1) - 2
            mask |= bit
        if valid and ell:
            result[ell] += 1
    return result


def bipartite_edges(left_size: int) -> tuple[tuple[int, int], ...]:
    return tuple((x, y) for x in range(left_size) for y in range(left_size, N))


def run_family(left_size: int) -> dict:
    pairs = bipartite_edges(left_size)
    total, nonzero, max_abs, total_starts, total_states, total_transitions = 1 << len(pairs), 0, 0, 0, 0, 0
    examples: list[dict] = []
    for edge_mask in range(total):
        edges = tuple(pair for j, pair in enumerate(pairs) if (edge_mask >> j) & 1)
        h, stats = h_by_subset_dp(adjacency(N, edges))
        r45 = h[3] - 3 * h[4] + 2 * h[5]
        total_starts += stats["ordered_nonedge_starts"]
        total_states += stats["nonzero_dp_cells"]
        total_transitions += stats["integer_transitions"]
        max_abs = max(max_abs, abs(r45))
        if r45:
            nonzero += 1
            if len(examples) < 8:
                examples.append({"edge_mask": edge_mask, "edges": [list(e) for e in edges],
                                 "h3": h[3], "h4": h[4], "h5": h[5], "R45": r45})
    return {"fixed_bipartition": [left_size, N - left_size], "edge_subsets": total,
            "nonzero_R45": nonzero, "max_abs_R45": max_abs, "counterexample_examples": examples,
            "dp_totals": {"ordered_nonedge_starts": total_starts, "nonzero_dp_cells": total_states,
                          "integer_transitions": total_transitions}}


def cross_checks() -> dict:
    # P7 exercises all target buckets and independently checks the rank-to-L
    # conversion against 7! direct words.
    p7_adj = adjacency(7, tuple((i, i + 1) for i in range(6)))
    p7_dp, _ = h_by_subset_dp(p7_adj)
    p7_words = h_by_direct_words(p7_adj)
    assert p7_dp == p7_words and (p7_dp[3], p7_dp[4], p7_dp[5]) == (128, 64, 32)
    # K_{3,4} is a nontrivial n=7 member of the enumerated family.
    k34 = tuple((x, y) for x in range(3) for y in range(3, 7))
    k34_adj = adjacency(7, k34)
    k34_dp, _ = h_by_subset_dp(k34_adj)
    k34_words = h_by_direct_words(k34_adj)
    assert k34_dp == k34_words
    # Boundary graph: K_{2,3} on {0,1}|{2,3,4}, with leaf 5 attached to 2
    # (the degree-two side; attaching it to the other side is a different graph).
    boundary_edges = tuple((x, y) for x in range(2) for y in range(2, 5)) + ((2, 5),)
    boundary_adj = adjacency(6, boundary_edges)
    b_dp, _ = h_by_subset_dp(boundary_adj)
    b_words = h_by_direct_words(boundary_adj)
    assert b_dp == b_words
    assert (b_dp[3], b_dp[4], b_dp[5]) == (16, 16, 0)
    assert b_dp[3] - 3 * b_dp[4] + 2 * b_dp[5] == -32
    return {"n7_P7_dp_equals_direct_word_reference": True,
            "n7_P7_h3_h4_h5": [p7_dp[3], p7_dp[4], p7_dp[5]],
            "n7_K34_dp_equals_direct_word_reference": True,
            "n7_K34_h3_h4_h5": [k34_dp[3], k34_dp[4], k34_dp[5]],
            "n6_K23_plus_leaf_edges": [list(e) for e in boundary_edges],
            "n6_K23_plus_leaf_dp_equals_direct_word_reference": True,
            "n6_K23_plus_leaf_h3_h4_h5": [b_dp[3], b_dp[4], b_dp[5]],
            "n6_K23_plus_leaf_R45": -32}


def main() -> None:
    started = datetime.now(timezone.utc)
    t0 = time.perf_counter()
    families = [run_family(k) for k in (1, 2, 3)]
    assert all(f["nonzero_R45"] == 0 and f["max_abs_R45"] == 0 for f in families)
    result = {"status": "passed", "N": N, "definition": "R45=h3-3h4+2h5", "method": "integer subset-DP over (mask,L), independent of global-normalizer permutation enumeration", "transition_rule": "from each nonadjacent ordered (r,p), append unused v iff v has a neighbor in mask; first p-neighbor has L=new complete rank minus 2", "families": families, "graphs_scanned_up_to_relabelling": sum(f["edge_subsets"] for f in families), "cross_checks": cross_checks(), "run": {"utc_started": started.isoformat(), "utc_finished": datetime.now(timezone.utc).isoformat(), "wall_seconds": time.perf_counter() - t0, "python": sys.version, "implementation": platform.python_implementation(), "platform": platform.platform(), "source_sha256": sha(Path(__file__))}}
    out = HERE / "r45-prefix-dp-n7.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
