"""Complete finite verification of the corrected second-peak insertion identity."""
from collections import Counter
from itertools import permutations
import json
import sys

from exact import grid, peaks, vertices


def main():
    n = 3
    xy, adj = grid(n, n)
    vertex_count = n * n
    full = (1 << vertex_count) - 1
    insertion = Counter()
    for v in range(vertex_count):
        remaining = full ^ (1 << v)
        neighbours = adj[v] & remaining
        remaining_vertices = list(vertices(remaining))
        for values in permutations(range(1, vertex_count)):
            labels = [0] * vertex_count
            for x, value in zip(remaining_vertices, values):
                labels[x] = value
            old_peaks = peaks(labels, adj, remaining)
            outside = [u for u in old_peaks if not (neighbours >> u & 1)]
            if len(outside) != 1:
                continue
            u = outside[0]
            max_neighbour = max(labels[x] for x in vertices(neighbours))
            insertion[(u, v)] += max(0, labels[u] - max_neighbour)

    direct = Counter()
    distance = Counter()
    for labels in permutations(range(1, vertex_count + 1)):
        peaks_now = peaks(labels, adj)
        if len(peaks_now) != 2:
            continue
        u = max(peaks_now, key=lambda x: labels[x])
        v = min(peaks_now, key=lambda x: labels[x])
        direct[(u, v)] += 1
        distance[abs(xy[u][0] - xy[v][0]) + abs(xy[u][1] - xy[v][1])] += 1

    assert direct == insertion
    output = {
        "model": "3x3 free-boundary nearest-neighbour grid",
        "identity_verified": "sum over ell on G-v with Peak(ell)\\N(v)={u} of max(ell(u)-M_v(ell),0)",
        "insertion_total": sum(insertion.values()),
        "direct_A2": sum(direct.values()),
        "distance_histogram": dict(sorted(distance.items())),
        "pairwise_match": True,
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
