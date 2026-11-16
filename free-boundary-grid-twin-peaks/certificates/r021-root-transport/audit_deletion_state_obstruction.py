#!/usr/bin/env python3
"""Exact P_4 check: the square candidate is not hereditary under deletion."""
from itertools import permutations
import json
from pathlib import Path


# This induced path lies in P_5 square.  Its first two vertices are a
# horizontal edge with x-coordinates 1 and 2.
H = [(1, 0), (2, 0), (2, 1), (3, 1)]
EDGES = {
    frozenset((H[0], H[1])),
    frozenset((H[1], H[2])),
    frozenset((H[2], H[3])),
}


def connected_prefix(word):
    seen = {word[0]}
    for v in word[1:]:
        if not any(frozenset((v, u)) in EDGES for u in seen):
            return False
        seen.add(v)
    return True


def main():
    counts = {v: 0 for v in H}
    for word in permutations(H):
        if connected_prefix(word):
            counts[word[0]] += 1
    assert counts == {H[0]: 1, H[1]: 3, H[2]: 3, H[3]: 1}
    rhs = ((5 - 1 - 1) / (1 + 1)) ** 2
    assert counts[H[1]] / counts[H[0]] > rhs
    payload = {
        "status": "passed",
        "ambient_square": 5,
        "induced_path": H,
        "root_counts": {str(v): counts[v] for v in H},
        "ratio_a_2_0_over_a_1_0": counts[H[1]] / counts[H[0]],
        "square_candidate_rhs_at_i_1": rhs,
    }
    out = Path(__file__).with_name("deletion-state-obstruction.json")
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
