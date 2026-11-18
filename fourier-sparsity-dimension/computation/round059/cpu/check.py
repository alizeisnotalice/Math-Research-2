"""Enumerate small positive affine covers using exact arithmetic."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import time


def span(vectors):
    space = {0}
    for vector in vectors:
        space |= {point ^ vector for point in space}
    return frozenset(space)


def affine_flats(n):
    size = 1 << n
    directions = {span(basis) for length in range(n + 1)
                  for basis in combinations(range(1, size), length)}
    flats = {}
    for direction in directions:
        normals = tuple(w for w in range(1, size)
                        if all((w & point).bit_count() % 2 == 0 for point in direction))
        for shift in range(size):
            points = tuple(sorted(shift ^ point for point in direction))
            flats[points] = normals
    return [(points, flats[points]) for points in sorted(flats)]


def binary_rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return len(pivots)


def modular_rank(rows, prime=65521):
    basis = {}
    for row in rows:
        row = [value % prime for value in row]
        for pivot, value in enumerate(row):
            if value == 0:
                continue
            if pivot not in basis:
                inverse = pow(value, -1, prime)
                basis[pivot] = [entry * inverse % prime for entry in row]
                break
            row = [(a - value * b) % prime for a, b in zip(row, basis[pivot])]
    return len(basis)


def census(n, target):
    size = 1 << n
    flats = affine_flats(n)
    by_point = [[index for index, (points, _) in enumerate(flats) if point in points]
                for point in range(size)]
    normals_by_flat = [normals for _, normals in flats]
    rows = [
        tuple((1 - 2 * ((w & points[0]).bit_count() % 2)) if w in normals else 0
              for w in range(1, size))
        for points, normals in flats
    ]
    profiles = Counter()
    seen = set()
    best = {}
    violations = {}
    witness = {}

    def record(assignment):
        support = tuple(sorted(assignment))
        if support in seen:
            return
        seen.add(support)
        assert all(sum(assignment[i] for i in by_point[point] if i in assignment) == target
                   for point in range(size))
        m = len(support)
        activity = [sum(i in assignment for i in by_point[point]) for point in range(size)]
        maximum = max(activity)
        t = binary_rank(w for index in support for w in normals_by_flat[index])
        r = modular_rank(rows[index] for index in support)
        profiles[(m, maximum, t, r)] += 1
        for name, gap in (("missed", m - maximum - t), ("centered", r - t)):
            if gap < best.get(name, size):
                best[name] = gap
            if gap < 0 and (name not in violations or m < violations[name]["m"]):
                violations[name] = {"m": m, "maximum": maximum, "t": t, "centered_rank": r,
                                    "activity": activity,
                                    "weighted_flats": [(flats[i][0], assignment[i]) for i in support]}
            if gap == 0 and (name not in witness or maximum > witness[name]["maximum"]
                             or (maximum == witness[name]["maximum"] and t > witness[name]["t"])):
                witness[name] = {"m": m, "maximum": maximum, "t": t, "centered_rank": r,
                                 "activity": activity,
                                 "weighted_flats": [(flats[i][0], assignment[i]) for i in support]}

    def complete(residual, assignment):
        if not any(residual):
            record(assignment)
            return
        point = next(i for i, value in enumerate(residual) if value)
        candidates = [index for index in by_point[point] if all(residual[i] for i in flats[index][0])]

        def choose(start):
            if residual[point] == 0:
                complete(residual, assignment)
                return
            for offset in range(start, len(candidates)):
                index = candidates[offset]
                domain = flats[index][0]
                limit = min(residual[i] for i in domain)
                for weight in range(1, limit + 1):
                    for i in domain:
                        residual[i] -= 1
                    assignment[index] = weight
                    choose(offset + 1)
                    del assignment[index]
                for i in domain:
                    residual[i] += limit

        choose(0)

    complete([target] * size, {})
    return {"n": n, "flat_catalog_size": len(flats), "integer_weight_target": target,
            "number_of_distinct_supports": len(seen), "minimum_gaps": best,
            "violations": violations, "equality_examples": witness,
            "profiles_m_maximum_t_centered_rank": {",".join(map(str, p)): count
                                                  for p, count in sorted(profiles.items())}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("target", type=int)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    result = census(args.n, args.target)
    result["elapsed_seconds"] = time.perf_counter() - started
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print({key: result[key] for key in ("n", "integer_weight_target", "number_of_distinct_supports",
                                      "minimum_gaps", "violations", "elapsed_seconds")})
