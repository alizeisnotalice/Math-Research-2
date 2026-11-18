"""Exact Boolean audit for overlapping four-relations and local rank two."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
from pathlib import Path
import platform
import time


def parity(value):
    return value.bit_count() & 1


def walsh(signed):
    signed = list(signed)
    stride = 1
    while stride < len(signed):
        for base in range(0, len(signed), 2 * stride):
            for point in range(base, base + stride):
                a, b = signed[point], signed[point + stride]
                signed[point], signed[point + stride] = a + b, a - b
        stride *= 2
    return signed


def rank(vectors):
    basis = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = vector
                break
            vector ^= basis[pivot]
    return len(basis)


def components(vertices, edges):
    origin = {vertex: vertex for vertex in vertices}

    def root(vertex):
        while origin[vertex] != vertex:
            vertex = origin[vertex]
        return vertex

    for edge in edges:
        representative = root(min(edge))
        for vertex in edge:
            origin[root(vertex)] = representative

    parts = {}
    for vertex in vertices:
        parts.setdefault(root(vertex), []).append(vertex)
    return list(parts.values())


def four_relations(labels):
    index = {label: i for i, label in enumerate(labels)}
    edges = []
    for a, b, c in combinations(sorted(labels), 3):
        d = a ^ b ^ c
        if d > c and d in index:
            edges.append([a, b, c, d])
    dimension = rank(sum(1 << index[x] for x in edge) for edge in edges)
    return edges, dimension


def selectors(a, shared):
    N = 1 << a
    data_dim = N + 1 if shared else 2 * N
    data_size = 1 << data_dim
    data_labels = [set() for _ in range(N)]
    table = []
    basis_planes = []
    for row in range(N):
        u = 1 if shared else 1 << (2 * row)
        v = 1 << (row + 1 if shared else 2 * row + 1)
        basis_planes.append((u, v))
    for data in range(data_size):
        for u, v in basis_planes:
            table.append(1 - 2 * (parity(data & u) * parity(data & v)))
    unique = {0}
    for row, (u, v) in enumerate(basis_planes):
        restricted = walsh(table[(point << a) | row] for point in range(data_size))
        assert {i: t for i, t in enumerate(restricted) if t} == {
            0: data_size // 2, u: data_size // 2, v: data_size // 2, u ^ v: -data_size // 2,
        }
        data_labels[row] = {i for i, value in enumerate(restricted) if i and value}
        unique.update(data_labels[row])
    expected = {(0, 0): len(table) // 2}
    for row, (u, v) in enumerate(basis_planes):
        for label, sign in ((u, 1), (v, 1), (u ^ v, -1)):
            for query in range(N):
                key = (query, label)
                expected[key] = expected.get(key, 0) + sign * (data_size // 2) * (1 - 2 * parity(query & row))
    actual = {(frequency & (N - 1), frequency >> a): value
              for frequency, value in enumerate(walsh(table)) if value}
    assert actual == {key: value for key, value in expected.items() if value}
    tail = len(actual) - 1
    expected_rank, expected_tail = (N + 1, 2 * N * N + 1) if shared else (2 * N, 3 * N * N)
    assert (rank(unique), tail) == (expected_rank, expected_tail)
    edges, u = four_relations(unique)
    assert u == N
    expected_edges = N * (N + 1) // 2 if shared else N
    assert len(edges) == expected_edges
    assert len(components(unique, edges)) == 1
    assert all(rank(row) == 2 for row in data_labels)
    return {"N": N, "shared_private_axis": shared, "input_dimension": a + data_dim,
            "truth_table_entries": len(table), "external_rank": rank(unique),
            "tail_points": tail, "four_relation_count": len(edges), "four_relation_rank": u,
            "four_relation_component_count": 1, "max_slice_linear_rank": 2,
            "tail_band_sizes": dict(sorted(Counter(label for _, label in actual if label).items())),
            "all_four_relations": edges}


def label_catalog(n):
    data_size = 1 << n
    catalog = []
    for mask in range(1 << data_size):
        spectrum = walsh(1 - 2 * ((mask >> point) & 1) for point in range(data_size))
        labels = tuple(i for i, value in enumerate(spectrum) if i and value)
        catalog.append((labels, rank(labels)))
    return catalog


def verify_chunk(n, begin, end, catalogs):
    ambient_size = 1 << n
    counter = Counter()
    for mask in range(begin, end):
        spectrum = walsh(1 - 2 * ((mask >> point) & 1) for point in range(ambient_size))
        for a in range(n + 1):
            query_size, data_size = 1 << a, 1 << (n - a)
            slices = []
            for u in range(query_size):
                local_mask = sum(((mask >> ((x << a) | u)) & 1) << x for x in range(data_size))
                labels, local_rank = catalogs[n - a][local_mask]
                if local_rank > 2:
                    break
                slices.append(labels)
            else:
                counter[(a, "accepted")] += 1
                bands = Counter(i >> a for i, value in enumerate(spectrum)
                                if i >> a and value)
                active = Counter(label for labels in slices for label in labels)
                assert bands.keys() == active.keys()
                assert all(bands[label] * count >= query_size for label, count in active.items())
                r = rank(bands)
                q = sum(bands.values())
                triangles = {labels for labels in slices if len(labels) == 3}
                assert all(len(labels) in (0, 1, 3) for labels in slices)
                assert all(triple[0] ^ triple[1] ^ triple[2] == 0 for triple in triangles)
                participating = {label for triple in triangles for label in triple}
                weights = {label: 2 if label in participating else 3 for label in bands}
                assert 3 * r <= sum(weights.values())
                assert all(sum(weights[label] ** 2 for label in labels) <= 12 for labels in slices)
                assert 3 * r * r <= 4 * q
                if triangles:
                    counter[(a, "has_triangles")] += 1
                if r and 3 * r * r == 4 * q:
                    counter[(a, "positive_rank_equality")] += 1
                    assert all(len(labels) == 3 for labels in slices)
                    assert len(set(active.values())) == 1
                    assert all(bands[label] * count == query_size for label, count in active.items())
                    parts = components(participating, triangles)
                    assert all(len(part) == 3 for part in parts)
                    assert r == 2 * len(parts)
    return counter


def exhaustive(workers):
    catalogs = [label_catalog(n) for n in range(5)]
    results = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for n in range(1, 5):
            size = 1 << (1 << n)
            work = [pool.submit(verify_chunk, n, size * i // workers, size * (i + 1) // workers, catalogs[:n+1])
                    for i in range(workers)]
            total = Counter()
            for item in work:
                total.update(item.result())
            results.append({"n": n, "truth_tables": size, "coordinate_splits": n + 1,
                            "by_query_dimension": [{k: total[(a,k)] for k in
                              ("accepted", "has_triangles", "positive_rank_equality")} for a in range(n + 1)]})
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    start = time.perf_counter()
    result = {"schema": "a28-r067-overlap-rank-two-audit-v1", "python": platform.python_version(),
              "workers": args.workers, "exact_small_census": exhaustive(args.workers),
              "boolean_selector_families": [selectors(a, shared) for a in range(3) for shared in (False, True)] + [selectors(3, True)],
              "elapsed_seconds": time.perf_counter() - start}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(result["exact_small_census"])
    print("runtime", result["elapsed_seconds"], "seconds")
