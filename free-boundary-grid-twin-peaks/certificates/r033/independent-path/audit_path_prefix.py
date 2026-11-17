#!/usr/bin/env python3
"""Independent audit of the first-merge path code on all P3-square words.

The analytic code is made explicit and decoded here.  The exhaustive n=3
check is only an injectivity audit; it is not used to infer the general bound.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from hashlib import sha256
from itertools import combinations, permutations
from math import comb
from pathlib import Path
import json
import platform
import sys
import time

HERE = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def grid_adjacency(side: int = 3) -> list[int]:
    n, adj = side * side, [0] * (side * side)
    for x in range(side):
        for y in range(side):
            u = x * side + y
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < side and 0 <= yy < side:
                    adj[u] |= 1 << (xx * side + yy)
    return adj


def one_peak(word: tuple[int, ...], adj: list[int]) -> bool:
    seen = 1 << word[0]
    for v in word[1:]:
        if not (adj[v] & seen):
            return False
        seen |= 1 << v
    return True


def double_peak(word: tuple[int, ...], adj: list[int]) -> bool:
    """Exactly two initial components: roots are nonadjacent, then connected growth."""
    p, r = word[:2]
    if adj[p] & (1 << r):
        return False
    seen = (1 << p) | (1 << r)
    for v in word[2:]:
        if not (adj[v] & seen):
            return False
        seen |= 1 << v
    return True


def parent_chain(v: int, root: int, parent: dict[int, int]) -> list[int]:
    chain = [v]
    while chain[-1] != root:
        chain.append(parent[chain[-1]])
    return chain


def encode(word: tuple[int, ...], adj: list[int]) -> tuple[tuple[int, ...], tuple] | None:
    """Return (one-peak output, explicit code), or None for non-double words.

    The parent of each pre-merge nonroot is its earliest already-exposed
    neighbour.  At the first merge, each branch neighbour is chosen by the
    same earliest-rank rule.  These choices make the code deterministic.
    """
    n = len(word)
    p, r = word[:2]
    if not double_peak(word, adj):
        return None
    rank = {v: i for i, v in enumerate(word)}
    seen = {p, r}
    branch = {p: 0, r: 1}
    parent: dict[int, int] = {}
    for pos in range(2, n):
        v = word[pos]
        neighbors = [u for u in seen if adj[v] & (1 << u)]
        if not neighbors:
            return None
        labels = {branch[u] for u in neighbors}
        if labels == {0, 1}:
            q, w = pos, v  # w has one-based rank q+1.
            a = min((u for u in neighbors if branch[u] == 0), key=rank.__getitem__)
            b = min((u for u in neighbors if branch[u] == 1), key=rank.__getitem__)
            left = list(reversed(parent_chain(a, p, parent)))
            right = parent_chain(b, r, parent)
            path = left + [w] + right  # oriented p -> w -> r
            assert len(path) == len(set(path)) and path[0] == p and path[-1] == r
            ell, j = len(path), path.index(w)
            assert 3 <= ell <= q + 1 and 1 <= j <= ell - 2
            path_set = set(path)
            output = tuple(path + [u for u in word if u not in path_set])
            if not one_peak(output, adj) or output[0] != p:
                seen_out = 1 << output[0]
                bad_position = None
                for out_pos, out_v in enumerate(output[1:], start=1):
                    if not (adj[out_v] & seen_out):
                        bad_position = out_pos
                        break
                    seen_out |= 1 << out_v
                raise AssertionError({"mapping_not_one_peak": True, "word": word, "q": q, "path": path,
                                      "output": output, "first_bad_output_zero_index": bad_position})
            assert set(output[:q + 1]) == set(word[:q + 1])
            slots = tuple(sorted(rank[u] for u in path if u not in {p, w}))
            left_slots = tuple(sorted(rank[u] for u in path[1:j]))
            assert len(slots) == ell - 2 and len(left_slots) == j - 1
            assert 1 in slots and 1 not in left_slots and all(1 <= x < q for x in slots)
            code = (q, ell, slots, j, left_slots)
            assert decode(output, code) == word
            return output, code
        # Before the first merge every legal extension attaches to exactly one
        # branch; the deterministic parent remains in that branch.
        assert len(labels) == 1
        branch[v] = next(iter(labels))
        parent[v] = min(neighbors, key=rank.__getitem__)
        seen.add(v)
    raise AssertionError("a connected graph word with two initial components must merge")


def decode(output: tuple[int, ...], code: tuple) -> tuple[int, ...]:
    """Inverse of the displayed code; it establishes the required injection."""
    q, ell, slots, j, left_slots = code
    n = len(output)
    path, p, w = output[:ell], output[0], output[j]
    assert len(slots) == ell - 2 and len(left_slots) == j - 1
    result: list[int | None] = [None] * n
    result[0], result[q] = p, w
    left_set = set(left_slots)
    for pos, vertex in zip(sorted(left_slots), path[1:j]):
        result[pos] = vertex
    right_slots = [pos for pos in slots if pos not in left_set]
    for pos, vertex in zip(right_slots, reversed(path[j + 1:])):
        result[pos] = vertex
    tail = iter(output[ell:])
    for pos, vertex in enumerate(result):
        if vertex is None:
            result[pos] = next(tail)
    try:
        next(tail)
        raise AssertionError("unconsumed output tail")
    except StopIteration:
        pass
    return tuple(result)  # type: ignore[arg-type]


def A2_capacity(q: int) -> int:
    return sum(comb(q - 1, ell - 2) * (2 ** (ell - 2) - 1) for ell in range(3, q + 2))


def B_capacity(q: int) -> int:
    return sum(comb(q - 2, ell - 3) * (2 ** (ell - 3)) for ell in range(3, q + 2))


def code_capacity_checks(n: int) -> dict[str, dict[str, int | bool]]:
    out = {}
    for q in range(2, n):
        a, b = A2_capacity(q), B_capacity(q)
        assert a == 3 ** (q - 1) - 2 ** (q - 1)
        assert b == 3 ** (q - 2)
        out[str(q)] = {"A2_code_space": a, "A2_formula": 3 ** (q - 1) - 2 ** (q - 1),
                       "B_code_space": b, "B_formula": 3 ** (q - 2), "identities_pass": True}
    return out


def main() -> None:
    started, t0 = datetime.now(timezone.utc), time.perf_counter()
    adj, n = grid_adjacency(), 9
    code_owner: dict[tuple[tuple[int, ...], tuple], tuple[int, ...]] = {}
    incoming: Counter[tuple[tuple[int, ...], int]] = Counter()
    per_q: Counter[int] = Counter()
    double_words = 0
    for word in permutations(range(n)):
        packed = encode(word, adj)
        if packed is None:
            continue
        double_words += 1
        output, code = packed
        q = code[0]
        key = (output, code)
        previous = code_owner.setdefault(key, word)
        assert previous == word, {"collision": True, "output": output, "code": code, "first": previous, "second": word}
        incoming[(output, q)] += 1
        per_q[q] += 1
    assert double_words == len(code_owner)
    maxima = {}
    for q in range(2, n):
        group_counts = [count for (output, qq), count in incoming.items() if qq == q]
        observed = max(group_counts, default=0)
        assert observed <= B_capacity(q) <= A2_capacity(q)
        maxima[str(q)] = {"maximum_observed_preimages_per_fixed_output": observed,
                           "B_capacity": B_capacity(q), "A2_capacity": A2_capacity(q), "passes": True}
    result = {"status": "passed", "graph": "P3 square (3x3 grid)", "vertices": n,
              "permutations_examined": 362880, "double_peak_words_encoded": double_words,
              "method": "deterministic earliest-parent first-merge path code plus exact decode",
              "injection": {"same_output_plus_code_collisions": 0, "decode_reconstructed_every_encoded_word": True,
                            "fixed_q_prefix_vertex_set_preserved_every_time": True},
              "analytic_code_capacity_checks": code_capacity_checks(n), "observed_fixed_output_multiplicity": maxima,
              "double_words_by_q": {str(q): per_q[q] for q in sorted(per_q)},
              "run": {"utc_started": started.isoformat(), "utc_finished": datetime.now(timezone.utc).isoformat(),
                      "wall_seconds": time.perf_counter() - t0, "python": sys.version,
                      "implementation": platform.python_implementation(), "platform": platform.platform(),
                      "source_sha256": digest(Path(__file__))}}
    out = HERE / "path-prefix-independent-r033.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
