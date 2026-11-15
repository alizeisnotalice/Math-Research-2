#!/usr/bin/env python3
"""Small independent exact checks of R63; no claim of universal 3.5 selection."""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
from pathlib import Path
import hashlib
import json
import platform
import sys
import time


def subset_data(points):
    n = len(points)
    out = {}
    for mask in range(1, 1 << n):
        rows = [points[i] for i in range(n) if mask >> i & 1]
        xs, ys = [p[0] for p in rows], [p[1] for p in rows]
        m = sum((p[2] for p in rows), F(0))
        box = min(xs), max(xs), min(ys), max(ys)
        d = max(box[1] - box[0], box[3] - box[2])
        out[mask] = (m, d * d, box)
    return out


def halo_rects(points, bits):
    scale = 1 << bits
    lower, upper = [], []
    seen = set()
    for mask, (m, d2, box) in subset_data(points).items():
        if d2 >= m:
            continue
        closure = sum(1 << i for i, (x, y, _) in enumerate(points)
                      if box[0] <= x <= box[1] and box[2] <= y <= box[3])
        if closure != mask or closure in seen:
            continue
        seen.add(closure)
        assert all((v * scale).denominator == 1 for v in box)
        a, b, c, d = [int(v * scale) for v in box]
        numerator = m.numerator * scale * scale
        denominator = 4 * m.denominator
        floor = isqrt(numerator // denominator)
        assert floor * floor * denominator <= numerator
        assert (floor + 1) * (floor + 1) * denominator > numerator
        for dest, radius in ((lower, max(0, floor - 1)), (upper, floor + 1)):
            r = (b - radius, a + radius, d - radius, c + radius)
            if r[0] < r[1] and r[2] < r[3]:
                dest.append(r)
    gamma = sum((r[1] - r[0]) * (r[3] - r[2]) for r in upper)
    gamma -= sum((r[1] - r[0]) * (r[3] - r[2]) for r in lower)
    return lower, upper, gamma


def union_area(rects):
    if not rects:
        return 0
    events = {}
    for i, (a, b, c, d) in enumerate(rects):
        events.setdefault(a, []).append((i, c, d, 1))
        events.setdefault(b, []).append((i, c, d, -1))
    active = {}
    area, last = 0, min(events)
    for x in sorted(events):
        intervals = sorted(active.values())
        length, end = 0, None
        for lo, hi in intervals:
            if end is None or lo > end:
                length += hi - lo
                end = hi
            elif hi > end:
                length += hi - end
                end = hi
        area += (x - last) * length
        for i, c, d, sign in events[x]:
            if sign == 1:
                active[i] = c, d
            else:
                del active[i]
        last = x
    assert not active
    return area


def frac(value):
    return str(value.numerator) + "/" + str(value.denominator)


def analyze(name, points, bits=32):
    data = subset_data(points)
    total = sum((p[2] for p in points), F(0))
    eligible = {s: item for s, item in data.items() if item[1] < item[0]}
    L = max(item[0] for item in eligible.values())
    masks = [s for s, item in eligible.items() if item[0] == L]
    square_cap_checks = 0
    for m, d2, _ in data.values():
        if d2 <= L:
            assert m <= L
            square_cap_checks += 1
    full_lo, full_hi, full_gamma = halo_rects(points, bits)
    area_scale = 1 << (2 * bits)
    deleted = []
    for mask in masks:
        _, _, (a, b, c, d) = data[mask]
        closure = sum(1 << i for i, (x, y, _) in enumerate(points)
                      if a <= x <= b and c <= y <= d)
        assert closure == mask
        residual = [p for i, p in enumerate(points) if not (mask >> i & 1)]
        collar_checks = 0
        for x, y, _ in residual:
            # Membership in the irrational closed collar is tested by squaring
            # its nonnegative one-sided coordinate distances.
            distance = max(b - x, x - a, d - y, y - c, F(0))
            assert distance * distance > L
            collar_checks += 1
        exterior_checks = 0
        for subset, (mass, _, _) in data.items():
            if subset & mask == 0:
                assert mass <= data[subset | mask][1] - L
                exterior_checks += 1
        width, height = b - a, d - c
        assert width * width < L and height * height < L
        simple_cheap = 16 * min(width * width, height * height) <= 9 * L
        if residual:
            rem_lo, rem_hi, rem_gamma = halo_rects(residual, bits)
        else:
            rem_lo, rem_hi, rem_gamma = [], [], 0
        low_int = union_area(full_lo + rem_hi) - union_area(rem_hi)
        high_int = union_area(full_hi + rem_lo) - union_area(rem_lo)
        assert 0 <= low_int <= high_int
        assert high_int - low_int <= full_gamma + rem_gamma
        low, high = F(low_int, area_scale), F(high_int, area_scale)
        local_indices = [i for i, (x, y, _) in enumerate(points)
                         if max(a-x, x-b, c-y, y-d, F(0)) ** 2 <= L]
        local = [points[i] for i in local_indices]
        local_residual = [points[i] for i in local_indices if not (mask >> i & 1)]
        assert max(m for m, d2, _ in subset_data(local).values() if d2 < m) == L
        assert sum((p[2] for p in local), F(0)) <= 9 * L
        loc_lo, loc_hi, _ = halo_rects(local, bits)
        lr_lo, lr_hi, _ = halo_rects(local_residual, bits)
        local_low = union_area(loc_lo + lr_hi) - union_area(lr_hi)
        local_high = union_area(loc_hi + lr_lo) - union_area(lr_lo)
        assert max(low_int, local_low) <= min(high_int, local_high)
        # Outward sqrt(L) bound for the elementary area envelope.
        sqrt_floor = isqrt((L.numerator << (2 * bits)) // L.denominator)
        sqrt_upper = F(sqrt_floor + 1, 1 << bits)
        envelope_upper = (width + sqrt_upper) * (height + sqrt_upper)
        assert low <= envelope_upper
        assert high <= envelope_upper + F(full_gamma + rem_gamma, area_scale)
        if simple_cheap:
            assert low < F(7, 2) * L
        deleted.append({
            "mask": mask, "proper": mask != (1 << len(points)) - 1,
            "mass": frac(L), "width": frac(width), "height": frac(height),
            "collar_outside_source_checks": collar_checks,
            "exterior_growth_subset_checks": exterior_checks,
            "localized_source_count": len(local),
            "localized_mass": frac(sum((p[2] for p in local), F(0))),
            "localization_delta_intervals_overlap": True,
            "thin_branch_proved": simple_cheap,
            "deletion_ratio_lower": frac(low / L),
            "deletion_ratio_upper": frac(high / L),
            "deletion_ratio_upper_display": float(high / L),
            "delta_band_width": frac(high - low),
            "gamma_sum": frac(F(full_gamma + rem_gamma, area_scale)),
            "certified_at_7over2": high <= F(7, 2) * L,
        })
    return {
        "name": name, "source_count": len(points),
        "points": [[frac(v) for v in p] for p in points],
        "total_mass": frac(total), "max_eligible_mass": frac(L),
        "all_nonempty_subsets_checked": len(data),
        "square_cap_checks": square_cap_checks,
        "maximizer_count": len(masks), "maximizers": deleted,
    }


def main():
    start = time.perf_counter()
    points = lambda rows: [tuple(F(v) for v in row) for row in rows]
    cases = [
        ("singleton", points([(0, 0, 1)])),
        ("separated", points([(0, 0, 1), (3, 0, 1)])),
        ("wide_chain", points([(i, 0, F(3, 4)) for i in range(5)])),
        ("two_clusters_bridge", points([
            (F(-17, 16), 0, F(3, 8)), (-1, F(1, 16), F(5, 8)),
            (0, 0, 1), (1, F(-1, 16), F(5, 8)), (F(17, 16), 0, F(3, 8))])),
    ]
    for spacing in (F(17, 16), F(19, 16), F(21, 16)):
        cases.append(("square_grid_" + frac(spacing).replace("/", "over"),
                      points([(spacing * i, spacing * j, F(1, 2))
                              for i in range(3) for j in range(3)])))
    cases.append(("asymmetric_frame", points([
        (0, 0, F(3, 8)), (1, 0, F(5, 8)), (2, 0, F(1, 2)),
        (0, 1, F(1, 2)), (2, 1, F(1, 2)),
        (0, 2, F(1, 2)), (1, 2, F(5, 8)), (2, 2, F(3, 8))])))
    cases.append(("bridge_with_irrelevant_far_sources", cases[3][1]
                  + points([(-10, 0, F(1, 4)), (10, 0, F(1, 4))])))
    result = {
        "status": "PASS exact finite checks; no universal 3.5 proof",
        "python": sys.version, "platform": platform.platform(),
        "arithmetic": "integer/Fraction; strict root brackets at 32 bits",
        "seed": None, "floating_decisions": False,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases": [analyze(name, p) for name, p in cases],
    }
    result["wall_seconds"] = time.perf_counter() - start
    output = Path(sys.argv[1])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "cases": len(result["cases"]),
                      "subsets": sum(c["all_nonempty_subsets_checked"] for c in result["cases"]),
                      "maximizers": sum(c["maximizer_count"] for c in result["cases"]),
                      "wall_seconds": result["wall_seconds"]}))


if __name__ == "__main__":
    main()
