#!/usr/bin/env python3
"""Build and exactly verify one connected nonterminal fixed-position R64 seed."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import hashlib
import json
import platform
import sys
import time


def text(q):
    q = F(q)
    return str(q.numerator) + "/" + str(q.denominator)


def main():
    started = time.perf_counter()
    alpha, delta, eta = F(15, 16), F(1, 64), F(1, 2048)
    a = alpha / 2
    near, far = 1 - a + delta, 1 + a - delta
    points, masses, labels = [], [], []
    for y in (-a, F(0), a):
        for x in (-a, F(0), a):
            points.append((x, y))
            masses.append(F(251, 1024) if x != 0 and y != 0 else F(1, 256))
            labels.append("S")
    for dx, dy in ((1, 0), (1, 1), (0, 1), (-1, 1),
                   (-1, 0), (-1, -1), (0, -1), (1, -1)):
        points.append((near * dx, near * dy))
        masses.append(F(1, 256))
        labels.append("near")
    for dx, dy in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        points.append((far * dx, far * dy))
        masses.append(F(3, 4))
        labels.append("far-corner")
    n, smask = len(points), (1 << 9) - 1
    assert sum(masses[:9]) == 1 and all(q > 0 for q in masses)
    xs, ys = sorted({p[0] for p in points}), sorted({p[1] for p in points})
    closures = {}
    for left, right in combinations_with_replacement(xs, 2):
        for bottom, top in combinations_with_replacement(ys, 2):
            mask = sum(1 << i for i, (x, y) in enumerate(points)
                       if left <= x <= right and bottom <= y <= top)
            if not mask or mask in closures:
                continue
            sub = [points[i] for i in range(n) if mask >> i & 1]
            width = max(p[0] for p in sub) - min(p[0] for p in sub)
            height = max(p[1] for p in sub) - min(p[1] for p in sub)
            d2 = max(width, height) ** 2
            mass = sum(masses[i] for i in range(n) if mask >> i & 1)
            bound = max(F(1), d2)
            assert mass <= bound, (mask, mass, bound)
            if d2 < 1 and mask != smask:
                assert mass <= 1 - eta, ("unique-S-margin", mask, mass)
            closures[mask] = (mass, d2, bound)
    eligible = [mask for mask, (mass, d2, _) in closures.items() if d2 < mass]
    maximal_mass = max(closures[mask][0] for mask in eligible)
    maximizers = [mask for mask in eligible if closures[mask][0] == maximal_mass]
    assert maximal_mass == 1 and maximizers == [smask]
    reached = 1
    while True:
        enlarged = reached
        for mask in eligible:
            if mask & enlarged:
                enlarged |= mask
        if enlarged == reached:
            break
        reached = enlarged
    assert reached == (1 << n) - 1
    connection_margin = F(1, 8192)
    connection_rows = []
    for i in range(9, n):
        candidates = [mask for mask in eligible if mask >> i & 1 and mask & smask]
        chosen = min(candidates, key=lambda mask: (mask.bit_count(), closures[mask][1], mask))
        mass, d2, _ = closures[chosen]
        assert mass >= d2 + connection_margin
        connection_rows.append({"outside_index": i, "mask": chosen,
                                "minimum_mass": text(d2 + connection_margin),
                                "bbox_span_squared": text(d2)})
    total = sum(masses)
    actual_span2 = (2 * far) ** 2
    assert actual_span2 > F(35, 33) * total
    assert actual_span2 - F(35, 33) * total > F(1, 16)
    source_box = (-a, a, -a, a)
    assert all(max(a-x, x+a, a-y, y+a, F(0)) > 1 for x, y in points[9:])
    assert all(max(-a-x, x-a, -a-y, y-a, F(0)) <= 1 for x, y in points)
    payload = {
        "schema": "r064-fixed-position-mass-polytope-v1",
        "status": "exactly verified connected nonterminal seed; no halo objective evaluated",
        "id": "alpha15over16-near8-farcorners4-positive21",
        "normalization": {"threshold": "1", "maximum_eligible_mass": "1"},
        "parameters": {"alpha": text(alpha), "delta": text(delta), "unique_margin_eta": text(eta)},
        "points": [[text(v) for v in p] for p in points],
        "initial_masses": [text(m) for m in masses],
        "source_roles": labels,
        "S_indices": list(range(9)),
        "S_mask": smask,
        "mass_constraints": {
            "nonnegative": True,
            "S_sum": "1/1",
            "full_bbox_caps": "sum(m_i for p_i in B) <= max(1, d(B)^2) for every source-coordinate bbox closure",
            "unique_max_design_restriction": "For each closure B different from S with d(B)^2<1, impose m(B)<=1-eta",
            "positive_floor_for_initial_chamber": "1/8192",
            "nonterminal_total_mass_cap": text(F(33, 35) * (actual_span2 - F(1, 16))),
            "connectivity_margin": text(connection_margin),
            "prescribed_connection_lower_constraints": connection_rows,
            "scope": "The floor, unique-max gap and nonterminal total cap are adversarial-design restrictions, not general theorems.",
        },
        "actual_seed_checks": {
            "positive_source_count": n,
            "unique_bbox_closures": len(closures),
            "strict_eligible_closures": len(eligible),
            "maximum_eligible_masks": maximizers,
            "eligible_connected": True,
            "actual_positive_bbox_span_squared": text(actual_span2),
            "total_mass": text(total),
            "nonterminal_margin_d2_minus_35M_over_33": text(actual_span2 - F(35, 33) * total),
            "all_sources_in_local_R": True,
            "outside_S_avoids_closed_collar": True,
        },
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integer/Fraction; no random numbers or floating decisions",
        "builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "wall_seconds": time.perf_counter() - started,
    }
    out = Path(sys.argv[1])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["actual_seed_checks"], indent=2))
    print("output_sha256", hashlib.sha256(out.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
