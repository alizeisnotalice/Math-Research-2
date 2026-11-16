#!/usr/bin/env python3
"""Direct 9! labeling enumeration for (distance, lower-peak birth rank)."""
from collections import Counter, defaultdict
from itertools import permutations
import json
from pathlib import Path


def main():
    n = 3
    N = n * n
    adj = []
    for y in range(n):
        for x in range(n):
            adj.append(tuple(yy*n+xx for xx, yy in ((x-1,y),(x+1,y),(x,y-1),(x,y+1))
                             if 0 <= xx < n and 0 <= yy < n))
    joint = defaultdict(Counter)
    distance = Counter()
    birth = Counter()
    for labels in permutations(range(N)):
        peaks = [v for v in range(N) if all(labels[v] > labels[u] for u in adj[v])]
        if len(peaks) != 2:
            continue
        q, p = sorted(peaks, key=labels.__getitem__, reverse=True)
        assert labels[q] == N-1
        tau = N - labels[p]
        assert tau >= 2
        d = abs(q % n - p % n) + abs(q // n - p // n)
        joint[d][tau] += 1
        distance[d] += 1
        birth[tau] += 1
    out = {
        "n": n, "N": N, "method": "all 9! label permutations",
        "joint_distance_birth": {str(d): {str(t): str(c) for t, c in sorted(row.items())}
                                 for d, row in sorted(joint.items())},
        "distance_counts": {str(d): str(c) for d, c in sorted(distance.items())},
        "birth_counts": {str(t): str(c) for t, c in sorted(birth.items())},
        "A2": str(sum(distance.values())),
    }
    path = Path(__file__).with_name("n3-direct.json")
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
