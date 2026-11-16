#!/usr/bin/env python3
"""Exact small-word audit for the strengthened early-interface lemma."""

from itertools import combinations, permutations
from math import comb, factorial


def binary_words(total, b):
    for places in combinations(range(total), b):
        positions = set(places)
        yield tuple("B" if i in positions else "H" for i in range(total))


def count_patterns(a, b, g):
    total = a + b + g
    good = bad = 0
    for word in binary_words(total, b):
        first_b = word.index("B")
        h_before = sum(letter == "H" for letter in word[:first_b])
        # H is A_1,...,A_a,C_1,...,C_g, so C_1 has H-index a.
        if h_before <= a:
            good += 1
        else:
            bad += 1
    assert bad == comb(b + g - 1, b), (a, b, g, bad)
    assert good == comb(total, b) - comb(b + g - 1, b)
    return good, bad


def neighbors(v):
    graph = {
        "a0": {"a1", "c1"},
        "a1": {"a0", "c0"},
        "b0": {"b1"},
        "b1": {"b0", "b2"},
        "b2": {"b1", "c0"},
        "c0": {"a1", "b2"},
        "c1": {"a0"},
    }
    return graph[v]


def births(word):
    seen, result = set(), 0
    for v in word:
        if not any(u in seen for u in neighbors(v)):
            result += 1
        seen.add(v)
    return result


def audit_births():
    # A and B are fixed one-peak words. H=A followed by arbitrary C order.
    a_word = ("a0", "a1")
    b_word = ("b0", "b1", "b2")
    c_word = ("c0", "c1")
    total, b = 7, 3
    checked = 0
    for pattern in binary_words(total, b):
        first_b = pattern.index("B")
        h_before = sum(letter == "H" for letter in pattern[:first_b])
        if h_before > len(a_word):
            continue
        ai = bi = 0
        word = []
        for letter in pattern:
            if letter == "B":
                word.append(b_word[bi])
                bi += 1
            else:
                if ai < len(a_word):
                    word.append(a_word[ai])
                else:
                    word.append(c_word[ai - len(a_word)])
                ai += 1
        assert births(word) == 2, (pattern, word)
        checked += 1
    assert checked == comb(total, b) - comb(3 + 2 - 1, 3)
    return checked


def grid_neighbors(v, n):
    x, y = v
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        q = (x + dx, y + dy)
        if 0 <= q[0] < n and 0 <= q[1] < n:
            yield q


def connected(vertices, n):
    start = next(iter(vertices))
    seen, stack = {start}, [start]
    while stack:
        v = stack.pop()
        for q in grid_neighbors(v, n):
            if q in vertices and q not in seen:
                seen.add(q)
                stack.append(q)
    return seen == vertices


def attachment_orders(seed, remaining, n):
    total = 0
    for order in permutations(remaining):
        seen = set(seed)
        legal = True
        for v in order:
            if not any(q in seen for q in grid_neighbors(v, n)):
                legal = False
                break
            seen.add(v)
        total += legal
    return total


def audit_full_grid_one_sided(n):
    vertices = {(x, y) for x in range(n) for y in range(n)}
    A = {(0, 0)}
    C = {(x, 0) for x in range(1, n - 1)} | {(0, 1)}
    B = vertices - A - C
    assert connected(A, n) and connected(B, n)
    assert connected(vertices - A, n) and connected(vertices - B, n)
    assert not any(q in B for v in A for q in grid_neighbors(v, n))
    bridges = {v for v in C if any(q in A for q in grid_neighbors(v, n))
               and any(q in B for q in grid_neighbors(v, n))}
    assert bridges == {(1, 0), (0, 1)}
    g = len(C)
    assert attachment_orders(A, C, n) == g
    assert attachment_orders(B, C, n) == factorial(g)
    suffixes = []
    for w in bridges:
        suffixes.append(attachment_orders(A | B | {w}, C - {w}, n))
    assert suffixes == [factorial(g - 1), factorial(g - 1)]
    return {"n": n, "g": g, "h_A": g, "h_B": factorial(g),
            "bridge_suffixes": suffixes}


if __name__ == "__main__":
    for values in ((1, 1, 1), (2, 3, 2), (4, 2, 3), (5, 4, 4)):
        print(values, count_patterns(*values))
    print({"birth_checked_patterns": audit_births()})
    for n in range(3, 9):
        print(audit_full_grid_one_sided(n))
