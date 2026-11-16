#!/usr/bin/env python3
"""Exact small check of the one-peak interface-root gluing inequality."""

from itertools import permutations
from math import comb, factorial
import json
from pathlib import Path
import sys
from collections import Counter, deque


def graph(width, height):
    adj = []
    for y in range(height):
        for x in range(width):
            ns = []
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < width and 0 <= yy < height:
                    ns.append(yy * width + xx)
            adj.append(tuple(ns))
    return adj


def onepeak_counts(width, height):
    adj = graph(width, height)
    root_counts = [0] * (width * height)
    for lab in permutations(range(width * height)):
        roots = [v for v, ns in enumerate(adj) if all(lab[v] > lab[u] for u in ns)]
        if len(roots) == 1:
            root_counts[roots[0]] += 1
    return root_counts


def is_onepeak(order, adj, root=0):
    """The order is high-to-low; a non-first vertex needs an earlier neighbor."""
    position = [0] * len(order)
    for i, vertex in enumerate(order):
        position[vertex] = i
    return order[0] == root and all(
        any(position[u] < position[v] for u in adj[v]) for v in order[1:]
    )


def shortest_path(adj, start, goal):
    predecessor = {start: None}
    queue = deque([start])
    while queue:
        v = queue.popleft()
        if v == goal:
            path = []
            while v is not None:
                path.append(v)
                v = predecessor[v]
            return tuple(reversed(path))
        for u in adj[v]:
            if u not in predecessor:
                predecessor[u] = v
                queue.append(u)
    raise AssertionError("disconnected graph")


def relocate_to_corner(order, adj, corner=0):
    path = shortest_path(adj, corner, order[0])
    path_set = set(path)
    return path + tuple(v for v in order if v not in path_set)


def independent_peak_family(n):
    interior = [(x, y) for x in range(1, n - 1) for y in range(1, n - 1)]
    classes = [[v for v in interior if (v[0] + 2 * v[1]) % 5 == r]
               for r in range(5)]
    family = max(classes, key=len)
    assert len(family) * 5 >= len(interior)
    assert all(abs(x - xx) + abs(y - yy) >= 3
               for i, (x, y) in enumerate(family)
               for xx, yy in family[i + 1:])
    return len(family)


def local_peaks(order, adj):
    position = [0] * len(order)
    for i, vertex in enumerate(order):
        position[vertex] = i
    return [v for v, ns in enumerate(adj) if all(position[v] < position[u] for u in ns)]


def interface_envelope_relocate(order, width, height):
    """Relocate the connected right-interface envelope to its top corner."""
    adj = graph(width, height)
    interface = {y * width + width - 1 for y in range(height)}
    outside = [v for v in local_peaks(order, adj) if v not in interface]
    assert len(outside) <= 1
    envelope = set(interface)
    q = outside[0] if outside else None
    if q is not None:
        x, y = q % width, q // width
        envelope.update(y * width + xx for xx in range(x, width))
    # Deterministic BFS is a connected order from the top-right corner.
    root = width - 1
    queue = deque([root])
    seen = {root}
    connected_order = []
    while queue:
        v = queue.popleft()
        connected_order.append(v)
        for u in adj[v]:
            if u in envelope and u not in seen:
                seen.add(u)
                queue.append(u)
    assert seen == envelope
    return tuple(connected_order + [v for v in order if v not in envelope]), q, len(envelope)


def main():
    left = onepeak_counts(2, 2)
    full = onepeak_counts(4, 2)
    # Interface roots: top-right corner of the left square, top-left of right.
    injection = comb(7, 3) * left[1] * left[0]
    assert injection <= sum(full)

    # Independently check deterministic shortest-path relocation on 3 x 2.
    relocation_adj = graph(3, 2)
    fibres = Counter()
    source_count = 0
    for order in permutations(range(6)):
        roots = [v for v, ns in enumerate(relocation_adj)
                 if all(order.index(v) < order.index(u) for u in ns)]
        if len(roots) == 1:
            source_count += 1
            image = relocate_to_corner(order, relocation_adj)
            assert is_onepeak(image, relocation_adj, root=0)
            fibres[image] += 1
    assert source_count == sum(onepeak_counts(3, 2))
    assert max(fibres.values()) <= 6 ** (3 + 2)
    independent_sizes = {str(n): independent_peak_family(n) for n in range(3, 17)}

    # Exact audit of the interface-envelope relocation on all 3 x 2 orders.
    envelope_adj = graph(3, 2)
    envelope_fibres = Counter()
    envelope_sources = 0
    for order in permutations(range(6)):
        interface = {2, 5}
        if sum(v not in interface for v in local_peaks(order, envelope_adj)) <= 1:
            envelope_sources += 1
            image, q, envelope_size = interface_envelope_relocate(order, 3, 2)
            assert is_onepeak(image, envelope_adj, root=2)
            envelope_fibres[image] += 1
            assert envelope_size <= 3 + 2 - 1
    assert max(envelope_fibres.values()) <= 7 * 6 ** (3 + 2 - 1)
    ans = {
        "two_by_two_A1": sum(left),
        "two_by_two_corner_root_counts": {"top_left": left[0], "top_right": left[1]},
        "four_by_two_A1": sum(full),
        "interface_gluing_injection_lower_count": injection,
        "normalized_factor_A_over_total": "1/2",
        "factorials": {"small": factorial(4), "large": factorial(8)},
        "status": "passed",
        "relocation_3_by_2": {
            "source_onepeak_count": source_count,
            "distinct_corner_images": len(fibres),
            "max_observed_fibre": max(fibres.values()),
            "proved_coarse_fibre_bound": 6 ** (3 + 2),
        },
        "independent_interior_peak_families_n_3_to_16": independent_sizes,
        "interface_envelope_3_by_2": {
            "sources_T_right": envelope_sources,
            "distinct_corner_images": len(envelope_fibres),
            "max_observed_fibre": max(envelope_fibres.values()),
            "proved_coarse_fibre_bound": 7 * 6 ** (3 + 2 - 1),
        },
    }
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("onepeak-gluing-audit.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps(ans, indent=2))


if __name__ == "__main__":
    main()
