#!/usr/bin/env python3
"""Verify a 4x4 witness where a reflected-root label swap loses a peak."""
import json

n = 4
# Row-major by y=0,1,2,3.  Larger labels are peaks when all neighbours are lower.
labels = [8, 2, 10, 5,
          14, 16, 15, 13,
          1, 6, 12, 3,
          9, 7, 11, 4]


def vertex(x, y):
    return y * n + x


def coord(v):
    return (v % n, v // n)


def neighbours(v):
    x, y = coord(v)
    return [vertex(xx, yy) for xx, yy in
            ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
            if 0 <= xx < n and 0 <= yy < n]


def peaks(a):
    return [v for v in range(n * n) if all(a[v] > a[u] for u in neighbours(v))]


p = vertex(0, 3)
r = vertex(3, 3)
after = labels.copy()
after[p], after[r] = after[r], after[p]
before_peaks = peaks(labels)
after_peaks = peaks(after)
assert sorted(map(coord, before_peaks)) == [(0, 3), (1, 1)]
assert [coord(v) for v in after_peaks] == [(1, 1)]
assert labels[p] == 9 and labels[vertex(2, 3)] == 11
assert after[r] == 9 < after[vertex(2, 3)]
print(json.dumps({
    'n': n,
    'labels_row_major': labels,
    'peaks_before': [coord(v) for v in before_peaks],
    'selected_peak': coord(p),
    'reflected_target': coord(r),
    'target_neighbour_with_larger_label': coord(vertex(2, 3)),
    'peaks_after_swap': [coord(v) for v in after_peaks],
}, indent=2))
