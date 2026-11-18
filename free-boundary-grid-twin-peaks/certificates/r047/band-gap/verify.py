"""Construct near inputs, move their second peaks, and check the full grid."""
from collections import deque
import json
import platform
import random


def neighbours(vertex, n):
    x, y = vertex
    return [
        (u, v)
        for u, v in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
        if 0 <= u < n and 0 <= v < n
    ]


def distance(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1])


def checked_peaks(labels, n):
    # Deliberately use independent edge enumeration for the peak audit.
    shadowed = set()
    for x in range(n):
        for y in range(n):
            a = (x, y)
            for b in ((x + 1, y), (x, y + 1)):
                if b in labels:
                    assert labels[a] != labels[b]
                    shadowed.add(a if labels[a] < labels[b] else b)
    return labels.keys() - shadowed


def near_word(path, q, target, n, rng):
    excluded = {q, *neighbours(q, n)}
    word = list(path)
    added = set(path)
    frontier = deque(path)
    assert not added & excluded
    while len(word) < target - 1:
        v = frontier.popleft()
        step = neighbours(v, n)
        rng.shuffle(step)
        for w in step:
            if w not in added and w not in excluded:
                word.append(w)
                added.add(w)
                frontier.append(w)
                if len(word) == target - 1:
                    break

    word.append(q)
    added.add(q)
    frontier = deque(word)
    while frontier:
        for w in neighbours(frontier.popleft(), n):
            if w not in added:
                word.append(w)
                added.add(w)
                frontier.append(w)
    assert len(word) == n * n
    return word


def check(n, r, rng):
    N = n * n
    b = r - 1
    J = 2 * r + 3
    p = (rng.randrange(n), rng.randrange(n))
    grid = {(x, y) for x in range(n) for y in range(n)}
    ball = {v for v in grid if distance(v, p) <= r}
    direction = 1 if p[0] < n // 2 else -1
    path = [(p[0] + i * direction, p[1]) for i in range(r + 2)]
    q = rng.choice(sorted(v for v in ball if distance(p, v) >= 3
                         and all(distance(w, v) >= 2 for w in path)))
    D = {(x, p[1]) for x in range(min(p[0], q[0]), max(p[0], q[0]) + 1)}
    D.update((q[0], y) for y in range(min(p[1], q[1]), max(p[1], q[1]) + 1))
    D.remove(q)
    S = {v for u in D for v in neighbours(u, n)} - D
    assert len(D) == distance(p, q) and len(S) < J and q in S
    target = rng.randint(N // 5, N // 2)
    word = near_word(path, q, target, n, rng)
    quantum = 1600 * J * (b + 1)
    ceiling = (N + 1) * quantum
    delta = ceiling // 800
    step = delta // J
    far_gap = ceiling // 200
    labels = {v: (N - i) * quantum for i, v in enumerate(word)}
    assert checked_peaks(labels, n) == {p, q}
    assert labels[p] > ceiling - delta and labels[q] < ceiling - 2 * delta
    assert sum(x > ceiling - far_gap for x in labels.values()) <= N // 100
    band = range(0, n // 4) if 2 * p[0] >= n else range(n - n // 4, n)
    candidates = [
        (x, y) for x in band for y in range(n) if (x + 2 * y) % 5 == 0
    ]

    verified = 0
    high_shell_cases = 0
    for j in range(J):
        lo = ceiling - 2 * delta + step * j
        hi = lo + step
        if any(lo <= labels[v] <= hi for v in S):
            continue
        H = {v for v in grid if labels[v] > hi}
        if H & S:
            high_shell_cases += 1
        assert p in H and q not in H
        assert all(
            any(labels[w] > labels[v] for w in neighbours(v, n))
            for v in H - {p}
        )

        B = D - H
        assert any(v in B for v in neighbours(q, n))
        seen = D & H
        queue = deque(sorted(seen))
        sequence = []
        while queue:
            for w in neighbours(queue.popleft(), n):
                if w in B and w not in seen:
                    seen.add(w)
                    sequence.append(w)
                    queue.append(w)
        assert len(sequence) == len(B)
        chosen = [
            (w, max(labels[v] for v in neighbours(w, n)))
            for w in candidates
            if labels[w] < hi
            and max(labels[v] for v in neighbours(w, n)) <= ceiling - far_gap
        ]
        assert len(chosen) * 40 >= N
        assert all(8 * distance(p, w) > n for w, _ in chosen)
        assert all(distance(p, v) > r + 1 for w in candidates for v in [w, *neighbours(w, n)])
        modified = labels.copy()
        for j, v in enumerate(sequence, 1):
            modified[v] = hi - j * (step // (len(B) + 1))
        assert all(lo < modified[v] < hi for v in B)
        assert checked_peaks(modified, n) == {p}
        w, maximum = rng.choice(chosen)
        modified[w] = (maximum + hi) // 2
        assert maximum < modified[w] < hi
        assert checked_peaks(modified, n) == {p, w}
        assert max(modified, key=modified.get) == p
        assert {v for v in grid if modified[v] > hi} == H
        verified += 1
    assert verified > 0
    return dict(n=n, radius=r, q_rank=target, empty_bands_checked=verified,
                high_shell_cases=high_shell_cases)


def main():
    seed = 47182026
    rng = random.Random(seed)
    cases = [check(n, r, rng) for n, r in [(64, 3), (80, 4), (96, 5), (128, 8)]]
    assert sum(item["high_shell_cases"] for item in cases) > 0
    print(json.dumps(dict(seed=seed, python=platform.python_version(), cases=cases), indent=2))


if __name__ == "__main__":
    main()
