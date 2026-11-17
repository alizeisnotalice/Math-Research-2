#!/usr/bin/env python3
"""Exact split of Psi into finite- and infinite-sole-child mechanisms."""
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
    assert len(w) == n * n and len(set(w)) == n * n
    return w


def saved_word():
    raw = SAVED.read_bytes()
    lines = raw.decode().splitlines()
    header = dict(zip(lines[0].split()[::2], lines[0].split()[1::2]))
    return int(header["n"]), list(map(int, lines[1].split()[1:])), hashlib.sha256(raw).hexdigest()


def split(n, w):
    N = n * n
    pos = [0] * N
    for i, v in enumerate(w, 1):
        pos[v] = i
    assert w[0] in range(N) and set(w) == set(range(N))
    assert all(any(pos[u] < pos[v] for u in nb(v, n)) for v in range(N) if v != w[0])
    finite = Fraction()
    infinity = Fraction()
    finite_count = infinity_count = 0
    for p in range(N):
        if p == w[0] or p in nb(w[0], n):
            continue
        ell = min(pos[u] for u in nb(p, n))
        b = min(
            (
                pos[v]
                for v in nb(p, n)
                if [u for u in nb(v, n) if pos[u] < pos[v]] == [p]
            ),
            default=N + 1,
        )
        R = b - ell - 1
        assert R > 0
        term = Fraction(ell - 1, R)
        if b == N + 1:
            infinity += term
            infinity_count += 1
        else:
            finite += term
            finite_count += 1
    return {
        "n": n,
        "finite_count": finite_count,
        "infinity_count": infinity_count,
        "psi_finite": [finite.numerator, finite.denominator],
        "psi_infinity": [infinity.numerator, infinity.denominator],
        "psi_total": [(finite + infinity).numerator, (finite + infinity).denominator],
    }


def add_comb_bound(row):
    n, N = row["n"], row["n"] * row["n"]
    T = (N + 2 * n - 1) // 2
    upper = Fraction(N * T, N - T)
    infinity = Fraction(*row["psi_infinity"])
    assert infinity <= upper
    row["comb_dominating_prefix_T"] = T
    row["psi_infinity_upper_bound"] = [upper.numerator, upper.denominator]
    row["psi_infinity_bound_pass"] = True
    return row


def main():
    n8, w8, sha = saved_word()
    out = {
        "kind": "exact-mixed-Psi-term-split",
        "python": platform.python_version(),
        "n8_saved_word_sha256": sha,
        "n8_low_search_word": split(n8, w8),
        "n9_row_word": split(9, row_word(9)),
        "comb_words": {
            str(n): add_comb_bound(split(n, comb_word(n))) for n in (5, 7, 9)
        },
        "scope": (
            "Deterministic words only. Exact integer/rational definition audit; "
            "no claim about all one-peak words or the one-peak law."
        ),
    }
    text = json.dumps(out, indent=2, sort_keys=True) + "\n"
    print(text, end="")
    print("output_sha256=" + hashlib.sha256(text.encode()).hexdigest())


if __name__ == "__main__":
    main()
