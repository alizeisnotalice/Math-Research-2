#!/usr/bin/env python3
"""Full 3-by-3 enumeration of four projection-skeleton classes.

This tests only the stronger proposed classwise injection into corner-root
orders.  A failure would not refute the original square root-ratio candidate.
"""
from collections import Counter, defaultdict
import json
from pathlib import Path


N = 3
VERTICES = [(x, y) for y in range(N) for x in range(N)]
ADJ = {
    v: tuple(
        w for w in VERTICES
        if abs(v[0] - w[0]) + abs(v[1] - w[1]) == 1
    )
    for v in VERTICES
}


def first_skeleton(word, coordinate):
    seen = {word[0][coordinate]}
    lo = hi = word[0][coordinate]
    ans = []
    for v in word[1:]:
        q = v[coordinate]
        if q in seen:
            continue
        assert q == lo - 1 or q == hi + 1
        if q == lo - 1:
            ans.append("L")
            lo = q
        else:
            ans.append("R")
            hi = q
        seen.add(q)
    assert len(ans) == N - 1
    return "".join(ans)


def last_skeleton(word, coordinate):
    counts = Counter(v[coordinate] for v in word)
    lo, hi = 0, N - 1
    ans = []
    for v in reversed(word):
        q = v[coordinate]
        counts[q] -= 1
        if counts[q]:
            continue
        if sum(counts.values()) == 0:
            continue
        assert q == lo or q == hi
        if q == lo:
            ans.append("L")
            lo += 1
        else:
            ans.append("R")
            hi -= 1
    assert len(ans) == N - 1
    return "".join(ans)


def skeleton(word):
    return (
        first_skeleton(word, 0),
        last_skeleton(word, 0),
        first_skeleton(word, 1),
        last_skeleton(word, 1),
    )


def words(root):
    out = []
    def walk(word, seen):
        if len(word) == N * N:
            out.append(tuple(word))
            return
        for v in VERTICES:
            if v not in seen and any(u in seen for u in ADJ[v]):
                walk(word + [v], seen | {v})
    walk([root], {root})
    return out


def main():
    root_words = {root: words(root) for root in VERTICES}
    corner = len(root_words[(0, 0)])
    assert corner == 1344
    report = {}
    maximum = (-1, None, None)
    for root, ws in root_words.items():
        classes = Counter(skeleton(w) for w in ws)
        for key, count in classes.items():
            if count > maximum[0]:
                maximum = (count, root, key)
        report[str(root)] = {
            "total_words": len(ws),
            "nonempty_classes": len(classes),
            "max_class_size": max(classes.values()),
            "classes_larger_than_corner": sum(c > corner for c in classes.values()),
        }
    payload = {
        "status": "passed",
        "n": N,
        "corner_root_count": corner,
        "maximum_class_size": maximum[0],
        "maximum_class_root": maximum[1],
        "maximum_class_skeleton": maximum[2],
        "all_classes_at_most_corner": maximum[0] <= corner,
        "by_root": report,
    }
    Path(__file__).with_name("n3-four-skeleton-audit.json").write_text(
        json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
