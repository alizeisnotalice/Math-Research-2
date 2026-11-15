"""Exact global p_d for n<=4 under uniform single-peak connected orders."""
from fractions import Fraction
import json
import sys

from exact import tables


def geodesic_prefixes(n, d):
    """All ordered d-vertex paths geodesic from their first vertex."""
    for root in range(n * n):
        x0, y0 = root % n, root // n

        def extend(path):
            if len(path) == d:
                yield path
                return
            vertex = path[-1]
            x, y = vertex % n, vertex // n
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                nx, ny = x + dx, y + dy
                if not (0 <= nx < n and 0 <= ny < n):
                    continue
                if abs(nx - x0) + abs(ny - y0) != len(path):
                    continue
                yield from extend(path + [ny * n + nx])

        yield from extend([root])


def result(n):
    _, _, one, _, completion, _ = tables(n, n, root_vectors=False)
    total = one[(1 << (n * n)) - 1]
    rows = []
    for d in range(1, 2 * n):
        sequences = list(geodesic_prefixes(n, d))
        numerator = sum(completion[sum(1 << vertex for vertex in path)] for path in sequences)
        p = Fraction(numerator, total)
        rows.append({
            "d": d,
            "geodesic_ordered_prefixes": len(sequences),
            "numerator": numerator,
            "p_d": str(p),
            "p_d_float": float(p),
            "geometric_path_prefix_mass": str(p * (2 ** (d - 2))) if d >= 2 else "1",
        })
    return {"n": n, "A1": total, "rows": rows}


if __name__ == "__main__":
    print(json.dumps({"definition": "global union over all roots and ordered graph-geodesic prefixes", "results": [result(2), result(3), result(4)]}, indent=2))
