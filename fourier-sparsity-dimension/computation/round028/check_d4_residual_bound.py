#!/usr/bin/env python3
"""Exact finite sanity checks for the R28 D4 column-count bound."""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
import random
import sys
import time
from pathlib import Path

OUT = Path(__file__).with_name("d4_residual_bound_result.json")


def rank_and_basis(words: list[int]) -> tuple[int, list[int]]:
    pivots: dict[int, int] = {}
    selected: list[int] = []
    for original in words:
        word = original
        while word:
            pivot = word.bit_length() - 1
            if pivot in pivots:
                word ^= pivots[pivot]
            else:
                pivots[pivot] = word
                selected.append(original)
                break
    return len(pivots), selected


def span(generators: list[int]) -> list[int]:
    code = [0]
    for generator in generators:
        code += [word ^ generator for word in code]
    return code


def column_statistics(weight4_basis: list[int], q: int) -> dict[str, object]:
    columns: dict[int, list[int]] = {}
    for point in range(q):
        label = sum(((row >> point) & 1) << i for i, row in enumerate(weight4_basis))
        columns.setdefault(label, []).append(point)
    d = len(weight4_basis)
    nonzero = {label: points for label, points in columns.items() if label}
    N = sum(len(points) for points in nonzero.values())
    r = sum(len(points) & 1 for points in nonzero.values())
    return {
        "d": d,
        "active_columns_N": N,
        "odd_nonzero_column_classes_r": r,
        "nonzero_class_sizes": sorted(len(points) for points in nonzero.values()),
        "N_bound_4d": 4 * d,
        "r_bound_floor_5d_over_2": (5 * d) // 2,
        "N_bound_passes": N <= 4 * d,
        "r_bound_passes": r <= (5 * d) // 2,
    }


def cubic_graph_model(name: str, vertices: list[int], edges: list[tuple[int, int]]) -> dict[str, object]:
    d, edge_count = len(vertices), len(edges)
    index = {vertex: i for i, vertex in enumerate(vertices)}
    edge_columns = [(1 << index[a]) | (1 << index[b]) for a, b in edges]
    matrix_columns = [1 << i for i in range(d)] + edge_columns
    row_weights = [sum((column >> i) & 1 for column in matrix_columns) for i in range(d)]
    # Embed the columns as distinct W-points, W=F_2^{edges}; no length-2^|E| word is built.
    edge_labels = [1 << j for j in range(edge_count)]
    vertex_labels = []
    for vertex in vertices:
        label = 0
        for j, edge in enumerate(edges):
            if vertex in edge:
                label ^= 1 << j
        vertex_labels.append(label)
    point_labels = vertex_labels + edge_labels
    row_label_xors = []
    for i in range(d):
        total = 0
        for column, label in zip(matrix_columns, point_labels):
            if (column >> i) & 1:
                total ^= label
        row_label_xors.append(total)
    rank, _ = rank_and_basis(matrix_columns)
    class_sizes: dict[int, int] = {}
    for column in matrix_columns:
        class_sizes[column] = class_sizes.get(column, 0) + 1
    return {
        "name": name,
        "d": d,
        "edges": edge_count,
        "ambient_W_dimension": edge_count,
        "ambient_q": 1 << edge_count,
        "row_weights": row_weights,
        "row_rank": rank,
        "matrix_column_count_N": len(matrix_columns),
        "all_matrix_columns_distinct": len(set(matrix_columns)) == len(matrix_columns),
        "all_physical_point_labels_distinct": len(set(point_labels)) == len(point_labels),
        "row_support_label_xors": row_label_xors,
        "RM1_orthogonality_verified": all(value == 0 for value in row_label_xors),
        "all_twin_classes_singleton": all(size == 1 for size in class_sizes.values()),
        "r": len(class_sizes),
        "five_d_over_two": 5 * d // 2,
    }


def random_affine_four_point_trials() -> dict[str, object]:
    rng = random.Random(20261117)
    q, trials, max_generators = 32, 80, 10
    outcomes: list[dict[str, int | bool]] = []
    for _ in range(trials):
        generators = []
        for _ in range(rng.randint(1, max_generators)):
            x = rng.randrange(q)
            a = rng.randrange(1, q)
            b = rng.randrange(1, q)
            while b == a:
                b = rng.randrange(1, q)
            plane = (1 << x) | (1 << (x ^ a)) | (1 << (x ^ b)) | (1 << (x ^ a ^ b))
            generators.append(plane)
        D_rank, independent_generators = rank_and_basis(generators)
        assert all(word.bit_count() == 4 for word in independent_generators)
        code = span(independent_generators)
        weight4 = [word for word in code if word.bit_count() == 4]
        d4, basis = rank_and_basis(weight4)
        assert all(word.bit_count() == 4 for word in basis)
        # The independent original generators already lie among these weight-four words.
        assert d4 == D_rank
        stats = column_statistics(basis, q)
        outcomes.append({
            "generator_count": len(generators),
            "D_rank": D_rank,
            "D_size": len(code),
            "D4_rank": d4,
            "N": int(stats["active_columns_N"]),
            "r": int(stats["odd_nonzero_column_classes_r"]),
            "N_ok": bool(stats["N_bound_passes"]),
            "r_ok": bool(stats["r_bound_passes"]),
        })
    return {
        "q": q,
        "trials": trials,
        "seed": 20261117,
        "max_input_generators": max_generators,
        "all_N_bounds_pass": all(item["N_ok"] for item in outcomes),
        "all_r_bounds_pass": all(item["r_ok"] for item in outcomes),
        "max_D4_rank": max(item["D4_rank"] for item in outcomes),
        "max_N": max(item["N"] for item in outcomes),
        "max_r": max(item["r"] for item in outcomes),
        "outcomes": outcomes,
    }


def main() -> None:
    started = time.perf_counter()
    k4_vertices = list(range(4))
    k4_edges = list(itertools.combinations(k4_vertices, 2))
    cube_vertices = list(range(8))
    cube_edges = [(v, v ^ (1 << bit)) for v in cube_vertices for bit in range(3) if v < (v ^ (1 << bit))]
    result = {
        "scope": "Exact finite checks only; the bound itself is proved combinatorially.",
        "python": sys.version,
        "platform": platform.platform(),
        "sharp_cubic_models": [
            cubic_graph_model("K4", k4_vertices, k4_edges),
            cubic_graph_model("cube_Q3", cube_vertices, cube_edges),
        ],
        "random_affine_weight4_trials": random_affine_four_point_trials(),
        "elapsed_seconds": time.perf_counter() - started,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if not all(model["RM1_orthogonality_verified"] and model["all_twin_classes_singleton"]
               for model in result["sharp_cubic_models"]):
        raise SystemExit(1)
    if not result["random_affine_weight4_trials"]["all_N_bounds_pass"] or not result["random_affine_weight4_trials"]["all_r_bounds_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
