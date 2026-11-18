from collections import deque
from fractions import Fraction
from itertools import product


def grid(width, height):
    points = frozenset(product(range(width), range(height)))
    neighbors = {
        v: points & {(v[0] - 1, v[1]), (v[0] + 1, v[1]), (v[0], v[1] - 1), (v[0], v[1] + 1)}
        for v in points
    }
    return points, neighbors


def complete(allowed, neighbors, prefix):
    used = set(prefix)
    queue = deque(prefix)
    ordering = list(prefix)
    while queue:
        v = queue.popleft()
        for next_vertex in sorted((neighbors[v] & allowed) - used):
            used.add(next_vertex)
            queue.append(next_vertex)
            ordering.append(next_vertex)
    assert used == allowed
    return ordering


def peaks_of_word(word, neighbors, allowed=None):
    peaks, seen = set(), set()
    for v in word:
        if allowed is None or v in allowed:
            if not (neighbors[v] & seen):
                peaks.add(v)
            seen.add(v)
    return peaks


def audit_forced_path():
    cases, splits = 0, 0
    for n in (5, 6, 9, 12, 20, 35):
        points, neighbors = grid(n, n)
        p, q = (0, 0), (n - 3, n - 3)
        path = [(x, 0) for x in range(n)] + [(n - 1, y) for y in range(1, n)]
        prefix = [p, q] + path[1:]
        rank = {v: k for k, v in enumerate(prefix)}
        for predecessor, v in zip(path, path[1:]):
            assert {a for a in neighbors[v] if a in rank and rank[a] < rank[v]} == {predecessor}
        assert not (neighbors[q] & set(path))
        word = complete(points, neighbors, prefix)
        assert peaks_of_word(word, neighbors) == {p, q}
        if n <= 6:
            for profile in product(range(1, n), repeat=n):
                for axis in (0, 1):
                    left = {v for v in points if v[axis] < profile[v[1 - axis]]}
                    assert len(peaks_of_word(word, neighbors, left)) != 1 or len(peaks_of_word(word, neighbors, points - left)) != 1
                    splits += 1
        cases += 1
    return cases, splits


def audit_reflection():
    records = []
    for n in (12, 14, 20, 26, 40, 60):
        points, neighbors = grid(n, n)
        half, local_neighbors = grid(n // 2, n)
        reflect = lambda v: (n - 1 - v[0], v[1])
        centers = set(product(range(2, n // 2 - 2, 5), range(2, n - 2, 5)))
        plus = {v for c in centers for v in local_neighbors[c] | {c}}
        core = half - plus
        base = complete(core, local_neighbors, [(0, 0)])
        left_high = {(x + 1, y) for x, y in centers}
        right_high = {(x - 1, y) for x, y in centers}
        word = (
            base + list(map(reflect, base))
            + sorted(left_high) + list(map(reflect, sorted(right_high)))
            + sorted(centers) + list(map(reflect, sorted(centers)))
            + sorted(plus - centers - left_high)
            + list(map(reflect, sorted(plus - centers - right_high)))
        )
        assert len(word) == len(set(word)) == len(points)
        assert peaks_of_word(word, neighbors) == {(0, 0), (n - 1, 0)}
        rank = {v: n * n - k for k, v in enumerate(word)}
        for max_on_left in (False, True):
            polarized = {}
            for v in half:
                w = reflect(v)
                hi, lo = max(rank[v], rank[w]), min(rank[v], rank[w])
                polarized[v], polarized[w] = (hi, lo) if max_on_left else (lo, hi)
            sorted_word = sorted(points, key=polarized.get, reverse=True)
            min_centers = {reflect(v) for v in centers} if max_on_left else centers
            assert peaks_of_word(sorted_word, neighbors) == min_centers | {(0, 0), (n - 1, 0)}
        records.append((n, len(centers) + 2))
    return records


def audit_bilateral_labels():
    records = []
    for n in (8, 12, 19, 27):
        points, neighbors = grid(n, n)
        split_x = n // 2
        left = {v for v in points if v[0] < split_x}
        right = points - left
        rows = tuple(range(3, n - 1, 5))
        left_new = {(split_x - 1, y) for y in rows}
        right_new = {(split_x, y) for y in rows}
        left_defect = {(split_x - 1, y - 1) for y in rows}
        right_defect = {(split_x, y - 1) for y in rows}
        left_low = ({(split_x - 2, y - 1) for y in rows} | {(split_x - 1, y - 2) for y in rows})
        right_low = ({(split_x + 1, y - 1) for y in rows} | {(split_x, y - 2) for y in rows})
        p, q = (0, 0), (n - 1, 0)
        left_high = complete(left - (left_new | left_defect | left_low), neighbors, [p])
        right_high = complete(right - (right_new | right_defect | right_low), neighbors, [q])
        word = left_high + right_high + sorted(left_new | right_new) + sorted(left_defect | right_defect) + sorted(left_low | right_low)
        assert len(word) == len(set(word)) == len(points)
        assert peaks_of_word(word, neighbors) == {p, q}
        assert peaks_of_word(word, neighbors, left - left_new) == {p} | left_defect
        assert peaks_of_word(word, neighbors, right - right_new) == {q} | right_defect
        assert peaks_of_word(word, neighbors, left) == {p}
        assert peaks_of_word(word, neighbors, right) == {q}
        for directions in product((False, True), repeat=len(rows)):
            profile = [split_x] * n
            for y, move_left in zip(rows, directions):
                profile[y] += 1 if move_left else -1
            side = {v for v in points if v[0] < profile[v[1]]}
            assert len(peaks_of_word(word, neighbors, side)) + len(peaks_of_word(word, neighbors, points - side)) >= 2 + len(rows)
        records.append((n, len(rows)))
    return records


def volumes_by_cell(ell, r, m_left, m_right):
    endpoints = sorted({Fraction(0), ell, r, m_left, m_right, Fraction(1)})
    intervals = list(zip(endpoints, endpoints[1:]))
    totals = [Fraction(0), Fraction(0), Fraction(0)]
    for i, (a, b) in enumerate(intervals):
        u = (a + b) / 2
        for j, (c, d) in enumerate(intervals):
            v = (c + d) / 2
            cell = (b - a) * (d - c)
            if m_left < u < ell and m_right < v < r:
                totals[0] += cell
            if not m_right and m_left < u < ell and i >= j:
                totals[1] += cell if i > j else cell / 2
            if not m_left and m_right < v < r and j >= i:
                totals[2] += cell if j > i else cell / 2
    return totals


def audit_volumes():
    cuts = [Fraction(i, 5) for i in range(5)]
    tested = 0
    for ell, r, m_left, m_right in product(cuts[1:], cuts[1:], cuts, cuts):
        independent = volumes_by_cell(ell, r, m_left, m_right)
        expected = [
            max(ell - m_left, 0) * max(r - m_right, 0),
            max(ell**2 - m_left**2, 0) * (not m_right) / 2,
            max(r**2 - m_right**2, 0) * (not m_left) / 2,
        ]
        assert independent == expected
        if expected[1] + expected[2] > 0:
            assert expected[0] > 0
        tested += 1
    return tested


if __name__ == "__main__":
    print("forced-path sizes and exhaustive axial profiles tested:", audit_forced_path())
    print("polarized n and exact peak count:", audit_reflection())
    print("bilateral defect n and active rows:", audit_bilateral_labels())
    print("exact two-dimensional volume cells, parameter tuples:", audit_volumes())
