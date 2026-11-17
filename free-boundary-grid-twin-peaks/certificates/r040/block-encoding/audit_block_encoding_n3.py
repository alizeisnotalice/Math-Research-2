#!/usr/bin/env python3
"""Full P3xP3 left-inverse audit for the low-branch block encoding."""
from collections import Counter
from itertools import permutations
import hashlib
import json
from math import comb
from pathlib import Path
import platform

HERE = Path(__file__).resolve().parent


def grid(n=3):
    adj = [set() for _ in range(n * n)]
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for xx, yy in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
                if 0 <= xx < n and 0 <= yy < n:
                    adj[v].add(yy * n + xx)
    return adj


def components(seen, adj):
    left = set(seen)
    ans = []
    while left:
        c = {left.pop()}
        todo = list(c)
        while todo:
            v = todo.pop()
            hit = adj[v] & left
            left -= hit
            c |= hit
            todo.extend(hit)
        ans.append(c)
    return ans


def two_peak_data(word, adj):
    seen = set()
    births = []
    merge = None
    for i, v in enumerate(word):
        cs = components(seen, adj) if seen else []
        hit = [c for c in cs if c & adj[v]]
        if not hit:
            births.append(v)
        if merge is None and len(hit) >= 2:
            merge = (i, v, hit)
        seen.add(v)
    if len(births) != 2:
        return None
    assert merge is not None
    q, w, hit = merge
    assert len(hit) == 2
    high = word[0]
    A = next(c for c in hit if high in c)
    B = next(c for c in hit if high not in c)
    assert len(A) + len(B) == q
    return q, w, A, B


def one_peak(word, adj):
    seen = set()
    for i, v in enumerate(word):
        if i and not (adj[v] & seen):
            return False
        seen.add(v)
    return True


def encode(word, adj):
    q, w, A, B = two_peak_data(word, adj)
    prefix = word[:q]
    Aorder = tuple(v for v in prefix if v in A)
    pi_B = tuple(v for v in prefix if v in B)
    J = tuple(i + 1 for i, v in enumerate(prefix) if v in B)
    position_B = {v: i + 1 for i, v in enumerate(pi_B)}
    x = min(B & adj[w])
    backward_path = [x]
    while backward_path[-1] != pi_B[0]:
        v = backward_path[-1]
        parent = min(
            (u for u in adj[v] & B if position_B[u] < position_B[v]),
            key=lambda u: position_B[u],
        )
        backward_path.append(parent)
    path_set = set(backward_path)
    P = tuple(sorted(position_B[v] for v in path_set))
    eta = tuple(backward_path) + tuple(v for v in pi_B if v not in path_set)
    assert P[0] == 1 and len(P) == len(backward_path)
    sigma = Aorder + (w,) + eta + word[q + 1 :]
    assert one_peak(sigma, adj)
    i = len(Aorder) + 1  # one-based position of w in sigma
    return sigma, i, P, J


def decode(code, b):
    sigma, i, P, J = code
    assert 2 <= i <= len(sigma) - b
    w = sigma[i - 1]
    Aorder = sigma[: i - 1]
    eta = sigma[i : i + b]
    suffix = sigma[i + b :]
    q = i + b - 1
    assert P and P[0] == 1 and all(1 <= j <= b for j in P)
    assert len(set(P)) == len(P)
    assert len(J) == b and all(2 <= j <= q for j in J)
    Jset = set(J)
    assert len(Jset) == len(J)
    path_forward = tuple(reversed(eta[: len(P)]))
    nonpath = iter(eta[len(P) :])
    path_by_pos = dict(zip(P, path_forward))
    pi_B = tuple(path_by_pos[j] if j in path_by_pos else next(nonpath) for j in range(1, b + 1))
    assert set(pi_B) == set(eta)
    assert tuple(nonpath) == ()
    ai = bi = 0
    prefix = []
    for j in range(1, q + 1):
        if j in Jset:
            prefix.append(pi_B[bi])
            bi += 1
        else:
            prefix.append(Aorder[ai])
            ai += 1
    assert ai == len(Aorder) and bi == len(pi_B)
    return tuple(prefix) + (w,) + tuple(suffix)


def path_adj(b):
    return [{u for u in (v - 1, v + 1) if 0 <= u < b} for v in range(b)]


def main():
    n, N = 3, 9
    adj = grid(n)
    A1 = A2 = 0
    low_counts = Counter()
    qb_counts = Counter()
    checked_codes = 0
    for word in permutations(range(N)):
        if one_peak(word, adj):
            A1 += 1
        data = two_peak_data(word, adj)
        if data is None:
            continue
        A2 += 1
        q, w, A, B = data
        assert word[0] in A and len(B) >= 1
        code = encode(word, adj)
        recovered = decode(code, len(B))
        assert recovered == word
        assert one_peak(code[0], adj)
        assert len(B) <= N - 2 and q <= N - 1
        low_counts[len(B)] += 1
        qb_counts[(q, len(B))] += 1
        checked_codes += 1
    assert A1 == 29568 and A2 == 143808
    rows = {}
    for b in range(1, N - 1):
        upper = A1 * (2 ** (b - 1)) * comb(N - 1, b + 1)
        rows[str(b)] = {
            "actual_low_branch_count": low_counts[b],
            "strong_encoding_upper_bound": upper,
            "bound_pass": low_counts[b] <= upper,
        }
    qb_rows = {}
    for q in range(2, N):
        for b in range(1, q):
            upper = A1 * (2 ** (b - 1)) * comb(q - 1, b)
            actual = qb_counts[(q, b)]
            assert actual <= upper
            qb_rows[f"q={q},b={b}"] = {
                "actual_count": actual,
                "q_b_encoding_upper_bound": upper,
                "bound_pass": True,
            }
    path_rows = {}
    for b in range(1, 10):
        adj_path = path_adj(b)
        all_onepeak = sum(one_peak(word, adj_path) for word in permutations(range(b)))
        fixed_endpoint = sum(
            one_peak(word, adj_path) for word in permutations(range(b)) if word[0] == 0
        )
        assert all_onepeak == 2 ** (b - 1) * fixed_endpoint and fixed_endpoint == 1
        path_rows[str(b)] = {
            "all_onepeak": all_onepeak,
            "fixed_endpoint_onepeak": fixed_endpoint,
            "factor_bound_equality": True,
        }
    result = {
        "kind": "full-P3xP3-low-branch-block-left-inverse",
        "python": platform.python_version(),
        "A1": A1,
        "A2": A2,
        "K2_codes_checked": checked_codes,
        "all_left_inverses_passed": checked_codes == A2,
        "rows_by_low_branch_size": rows,
        "rows_by_first_merge_size_and_low_branch": qb_rows,
        "path_fixed_root_sharpness_b_1_to_9": path_rows,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(text, end="")
    print("output_sha256=" + hashlib.sha256(text.encode()).hexdigest())


if __name__ == "__main__":
    main()
