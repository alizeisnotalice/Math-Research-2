#!/usr/bin/env python3
"""Exact audits for the boundary--sole-child bridge on specified one-peak words."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform

HERE = Path(__file__).resolve().parent
SAVED = HERE.parents[1] / 'r036/finite-candidates/n8-bfs-center-s2.txt'


def nb(v, n):
    x, y = v % n, v // n
    return [
        yy * n + xx
        for xx, yy in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
        if 0 <= xx < n and 0 <= yy < n
    ]


def row_word(n):
    return list(range(n * n))


def comb_word(n):
    assert n >= 5 and n % 2 == 1
    w = [y * n for y in range(n)]
    for y in range(0, n, 2):
        w.extend(y * n + x for x in range(1, n))
    w.extend(y * n + x for y in range(1, n, 2) for x in range(1, n))
    return w


def saved_word():
    raw = SAVED.read_bytes()
    lines = raw.decode().splitlines()
    h = dict(zip(lines[0].split()[::2], lines[0].split()[1::2]))
    return int(h["n"]), list(map(int, lines[1].split()[1:])), hashlib.sha256(raw).hexdigest()


def audit(n, w):
    N = n * n
    pos = [0] * N
    for i, v in enumerate(w, 1):
        pos[v] = i
    root = w[0]
    degree_root = len(nb(root, n))
    assert set(w) == set(range(N))
    early = {
        v: [u for u in nb(v, n) if pos[u] < pos[v]]
        for v in range(N)
    }
    assert all(early[v] for v in range(N) if v != root)
    eligible = {p for p in range(N) if p != root and p not in nb(root, n)}
    ell, b, term = {}, {}, {}
    for p in eligible:
        ell[p] = min(pos[u] for u in nb(p, n))
        b[p] = min(
            (pos[v] for v in nb(p, n) if early[v] == [p]),
            default=N + 1,
        )
        R = b[p] - ell[p] - 1
        assert R > 0
        term[p] = Fraction(ell[p] - 1, R)
    seen = set()
    singleton_children = 0
    checked = []
    for t, v in enumerate(w, 1):
        seen.add(v)
        if v != root and len(early[v]) == 1:
            singleton_children += 1
        boundary = sum(1 for u in seen for z in nb(u, n) if z not in seen)
        closed = seen | {z for u in seen for z in nb(u, n)}
        D = N - len(closed)
        finite_parents = {p for p in eligible if b[p] <= t}
        m = len(finite_parents)
        s_lower_numerator = N - t - D - degree_root
        parent_lower = max(0, (singleton_children + 3) // 4 - (1 + degree_root))
        weighted_parent_sum = sum((term[p] for p in finite_parents), Fraction())
        weighted_lower = Fraction(m * (m + 4), 8 * t)
        assert boundary <= degree_root + 2 * singleton_children
        assert boundary >= N - t - D
        assert singleton_children * 2 >= s_lower_numerator
        assert m >= parent_lower
        assert weighted_parent_sum >= weighted_lower
        if t in {1, N // 4, N // 2, 3 * N // 4, N}:
            checked.append(
                {
                    "t": t, "B_t": boundary, "D_t": D, "s_t": singleton_children,
                    "finite_eligible_parents": m,
                    "finite_weight_sum": [
                        weighted_parent_sum.numerator, weighted_parent_sum.denominator
                    ],
                    "finite_weight_lower": [
                        weighted_lower.numerator, weighted_lower.denominator
                    ],
                }
            )
    return {
        "n": n, "N": N, "root": root, "degree_root": degree_root,
        "eligible_count": len(eligible), "checked_times": checked,
        "all_prefix_bridge_assertions_passed": True,
    }


def main():
    n8, w8, sha = saved_word()
    out = {
        "kind": "exact-prefix-boundary-solechild-bridge",
        "python": platform.python_version(),
        "n8_saved_word_sha256": sha,
        "n8_low_search_word": audit(n8, w8),
        "n9_row_word": audit(9, row_word(9)),
        "n9_comb_word": audit(9, comb_word(9)),
        "scope": "Specified deterministic words only; no general asymptotic assertion.",
    }
    text = json.dumps(out, indent=2, sort_keys=True) + "\n"
    print(text, end="")
    print("output_sha256=" + hashlib.sha256(text.encode()).hexdigest())


if __name__ == "__main__":
    main()
