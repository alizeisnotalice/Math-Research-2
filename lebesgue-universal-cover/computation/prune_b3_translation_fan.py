#!/usr/bin/env python3
"""Exact coarse pruning in the B3 translation plane.

This is intentionally a *projection* of the 20-dimensional F8 placement
problem.  It uses only the fixed radius-1/2 disk and a rational triangle
strictly inside the normalized Reuleaux triangle B3.  Thus, when a cell is
pruned, every F8 placement whose B3 translation is in that cell has hull
area > tau.  A retained cell says nothing: the other 18 parameters have not
been searched or optimized.

All computational predicates below are integer predicates after scaling by
the denominator used for the disk polygon in outer_domain_certificate.json.
There is no floating-point arithmetic in either the pruning or the B3
membership checks.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from collections import Counter, deque
from pathlib import Path
from fractions import Fraction as F


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path(__file__).with_name("outer_domain_certificate.json")
TAU = F(419, 500)


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def tagged_hull(points):
    """CCW strict convex hull; each item is (x, y, material-label)."""
    # Coordinates are unique for the cells examined here.  Retaining the
    # first label is harmless if a rare equality occurs.
    by_xy = {}
    for x, y, tag in points:
        by_xy.setdefault((x, y), tag)
    ordered = [(x, y, by_xy[(x, y)]) for x, y in sorted(by_xy)]
    if len(ordered) < 3:
        raise AssertionError("not enough distinct points")

    def build(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and cross(h[-2], h[-1], p) <= 0:
                h.pop()
            h.append(p)
        return h[:-1]

    return build(ordered) + build(reversed(ordered))


def max_abs_interval(lo, hi):
    return max(abs(lo), abs(hi))


def check_inner_triangle():
    """Certify Q is in B3 = intersection of its three unit disks.

    The arc centres for B3 in the project convention are
       (1/sqrt(3),0), (-1/(2sqrt(3)), +/-1/2).
    We use a rational enclosure for sqrt(3) and interval upper bounds for
    the three squared distances.  The coordinate scale is imposed later.
    """
    q = [(F(11, 20), F(0)), (F(-11, 40), F(47, 100)), (F(-11, 40), F(-47, 100))]
    s_lo = F(1732050807568877, 10**15)
    s_hi = F(1732050807568878, 10**15)
    assert s_lo > 0 and s_lo * s_lo < 3 < s_hi * s_hi
    inv_lo, inv_hi = 1 / s_hi, 1 / s_lo
    centres = [
        ((inv_lo, inv_hi), (F(0), F(0))),
        ((-inv_hi / 2, -inv_lo / 2), (F(1, 2), F(1, 2))),
        ((-inv_hi / 2, -inv_lo / 2), (F(-1, 2), F(-1, 2))),
    ]
    upper = []
    for x, y in q:
        row = []
        for (cxlo, cxhi), (cylo, cyhi) in centres:
            dx = max_abs_interval(x - cxhi, x - cxlo)
            dy = max_abs_interval(y - cyhi, y - cylo)
            value = dx * dx + dy * dy
            assert value < 1
            row.append(value)
        upper.append(row)
    return q, {"sqrt3_interval": [str(s_lo), str(s_hi)], "squared_distance_upper_bounds": [[str(z) for z in r] for r in upper]}


def cell_record(cell, status, reason, min_area2=None, hull_vertices=None):
    level, x0, x1, y0, y1 = cell
    out = {
        "level": level,
        "x_scaled": [str(x0), str(x1)],
        "y_scaled": [str(y0), str(y1)],
        "status": status,
        "reason": reason,
    }
    if min_area2 is not None:
        out["minimum_twice_fan_area_scaled"] = str(min_area2)
    if hull_vertices is not None:
        out["fan_vertices_at_cell_center"] = hull_vertices
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--initial-level", type=int, default=3,
                        help="begin with 2^level by 2^level cells")
    parser.add_argument("--max-level", type=int, default=7)
    parser.add_argument("--node-budget", type=int, default=25000)
    args = parser.parse_args()
    if not (0 <= args.initial_level <= args.max_level <= 12):
        raise SystemExit("require 0 <= initial-level <= max-level <= 12")

    began = time.monotonic()
    raw = args.input.read_bytes()
    inp = json.loads(raw)
    den = int(inp["disk_polygon_denominator"])
    disk = [(int(x), int(y), "D") for x, y in inp["disk_polygon_points"]]
    assert all(x * x + y * y <= den * den // 4 for x, y, _ in disk)
    # At max level 7 translations have denominator 100*2^7, a divisor of
    # 10^15.  This keeps all coordinates as one common integer scale.
    assert den % (100 * (2 ** args.max_level)) == 0
    q, q_certificate = check_inner_triangle()
    q_scaled = [(int(x * den), int(y * den), "Q0" if i == 0 else f"Q{i}")
                for i, (x, y) in enumerate(q)]
    assert all(F(x, den) == q[i][0] and F(y, den) == q[i][1]
               for i, (x, y, _) in enumerate(q_scaled))

    # t_B3 in [-21/100,21/100]^2, from the verified round-14 global box.
    lo, hi = -21 * den // 100, 21 * den // 100
    assert hi - lo == 42 * den // 100
    threshold2 = 2 * TAU * den * den
    assert threshold2.denominator == 1
    threshold2 = threshold2.numerator

    q_by_tag = {tag: (x, y) for x, y, tag in q_scaled}

    def material_point(v, tx, ty):
        x, y, tag = v
        # The hull is formed at the cell centre, so a Q vertex stored in h
        # already carries that centre translation.  Recover its fixed Q base
        # by tag before applying the corner translation.  Using x+tx here
        # would incorrectly test Q+centre+t.
        return (q_by_tag[tag][0] + tx, q_by_tag[tag][1] + ty) if tag != "D" else (x, y)

    def assess(cell):
        level, x0, x1, y0, y1 = cell
        assert (x0 + x1) % 2 == 0 and (y0 + y1) % 2 == 0
        center = ((x0 + x1) // 2, (y0 + y1) // 2)
        h = tagged_hull(disk + [(x + center[0], y + center[1], tag) for x, y, tag in q_scaled])
        # Guard against a second application of the centre translation.
        assert all(material_point(v, *center) == (v[0], v[1]) for v in h)
        corners = [(x0, y0), (x0, y1), (x1, y0), (x1, y1)]
        areas = []
        for tx, ty in corners:
            vertices = [material_point(v, tx, ty) for v in h]
            # Strict positivity at four corners makes every fan determinant
            # positive throughout the rectangle because each is affine in t.
            if any(det(a, b) <= 0 for a, b in zip(vertices, vertices[1:] + vertices[:1])):
                return "retain", "fan_orientation_not_certified", None, len(h)
            areas.append(sum(det(a, b) for a, b in zip(vertices, vertices[1:] + vertices[:1])))
        # Every summand (and hence the shoelace sum) is affine: the possible
        # t,t determinant is zero.  Check this identity explicitly instead
        # of relying on that algebra in the program logic.
        if areas[0] + areas[3] != areas[1] + areas[2]:
            return "retain", "fan_area_affinity_check_failed", None, len(h)
        amin = min(areas)
        if amin > threshold2:
            return "prune", "strict_inner_fan_area_exceeds_tau", amin, len(h)
        return "retain", "inner_fan_lower_bound_not_above_tau", amin, len(h)

    queue = deque()
    n0 = 2 ** args.initial_level
    step = (hi - lo) // n0
    assert step * n0 == hi - lo and step % 2 == 0
    for i in range(n0):
        for j in range(n0):
            queue.append((args.initial_level, lo + i * step, lo + (i + 1) * step,
                          lo + j * step, lo + (j + 1) * step))

    leaves = []
    processed = 0
    counts = Counter()
    while queue and processed < args.node_budget:
        cell = queue.popleft()
        processed += 1
        status, reason, amin, hlen = assess(cell)
        counts[(status, reason)] += 1
        level, x0, x1, y0, y1 = cell
        if status == "prune":
            leaves.append(cell_record(cell, status, reason, amin, hlen))
        elif level < args.max_level:
            mx, my = (x0 + x1) // 2, (y0 + y1) // 2
            queue.extend(((level + 1, x0, mx, y0, my), (level + 1, x0, mx, my, y1),
                          (level + 1, mx, x1, y0, my), (level + 1, mx, x1, my, y1)))
        else:
            leaves.append(cell_record(cell, status, reason, amin, hlen))

    # A budget stop is explicitly retained rather than silently discarded.
    while queue:
        cell = queue.popleft()
        counts[("retain", "node_budget_not_processed")] += 1
        leaves.append(cell_record(cell, "retain", "node_budget_not_processed"))

    pruned = [r for r in leaves if r["status"] == "prune"]
    retained = [r for r in leaves if r["status"] == "retain"]
    scale_area = den * den
    # Each terminal rectangle has area (dx dy)/den^2.  This is a statement
    # about the B3 translation plane only, not a 20-dimensional volume.
    def plane_area(records):
        return sum((int(r["x_scaled"][1]) - int(r["x_scaled"][0])) *
                   (int(r["y_scaled"][1]) - int(r["y_scaled"][0])) for r in records)

    result = {
        "status": "completed_exact_coarse_projection_pruning",
        "scope": (
            "Rigorous only for the two B3 translation coordinates after the round-14 "
            "normalization. A pruned rectangle cannot occur in any area <419/500 F8 "
            "placement because conv(D union B3) already has an inner fan of area >tau. "
            "Retained rectangles are not feasible witnesses and no full 20D branch-and-bound "
            "or universal-cover lower bound is claimed."
        ),
        "source_global_box": {"input": args.input.name, "sha256": hashlib.sha256(raw).hexdigest()},
        "tau": str(TAU),
        "coordinate_denominator": str(den),
        "translation_domain_scaled": [str(lo), str(hi)],
        "translation_domain": "[-21/100,21/100]^2",
        "initial_level": args.initial_level,
        "max_level": args.max_level,
        "node_budget": args.node_budget,
        "processed_nodes": processed,
        "budget_exhausted": processed >= args.node_budget,
        "inner_material": {
            "disk": "the 512 rational points certified inside D in outer_domain_certificate.json",
            "B3_triangle": [[str(x), str(y)] for x, y in q],
            "B3_membership_certificate": q_certificate,
        },
        "criterion": {
            "fan": "center-cell convex hull of disk points union (Q+t), transported to cell corners",
            "orientation": "each consecutive determinant is affine in t and is strictly positive at all four corners",
            "area": "twice fan area is verified affine by A00+A11=A01+A10; its rectangular minimum is a corner",
            "prune_condition": "minimum twice fan area > 2*tau*denominator^2",
        },
        "counts_by_assessment": {f"{a}:{b}": n for (a, b), n in sorted(counts.items())},
        "terminal_leaf_counts": {"pruned": len(pruned), "retained": len(retained), "all": len(leaves)},
        "translation_plane_area": {
            "meaning": "areas in the original B3 translation square [-21/100,21/100]^2; fractions below are relative to that square, not absolute whole-plane areas",
            "domain_scaled_numerator": str((hi - lo) ** 2),
            "pruned_scaled_numerator": str(plane_area(pruned)),
            "retained_scaled_numerator": str(plane_area(retained)),
            "common_denominator": str(scale_area),
            "fraction_of_original_translation_square": {
                "pruned": str(F(plane_area(pruned), (hi - lo) ** 2)),
                "retained": str(F(plane_area(retained), (hi - lo) ** 2)),
            },
        },
        "leaves": leaves,
        "python": sys.version,
        "platform": platform.platform(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.monotonic() - began,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output), "processed_nodes": processed,
        "budget_exhausted": result["budget_exhausted"],
        "terminal_leaf_counts": result["terminal_leaf_counts"],
        "translation_plane_area": result["translation_plane_area"],
        "elapsed_seconds": result["elapsed_seconds"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
