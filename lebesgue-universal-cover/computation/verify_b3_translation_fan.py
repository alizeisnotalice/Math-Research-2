#!/usr/bin/env python3
"""Independent exact verifier for prune_b3_translation_fan.py output."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path


TAU = F(419, 500)


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def hull(points):
    labels = {}
    for x, y, label in points:
        labels.setdefault((x, y), label)
    p = [(x, y, labels[(x, y)]) for x, y in sorted(labels)]
    def half(seq):
        out = []
        for q in seq:
            while len(out) >= 2 and cross(out[-2], out[-1], q) <= 0:
                out.pop()
            out.append(q)
        return out[:-1]
    return half(p) + half(reversed(p))


def inner_triangle():
    q = [(F(11, 20), F(0)), (F(-11, 40), F(47, 100)), (F(-11, 40), F(-47, 100))]
    lo, hi = F(1732050807568877, 10**15), F(1732050807568878, 10**15)
    assert lo * lo < 3 < hi * hi
    ilo, ihi = 1 / hi, 1 / lo
    centres = [((ilo, ihi), (F(0), F(0))), ((-ihi / 2, -ilo / 2), (F(1, 2), F(1, 2))),
               ((-ihi / 2, -ilo / 2), (F(-1, 2), F(-1, 2)))]
    for x, y in q:
        for (xlo, xhi), (ylo, yhi) in centres:
            dx, dy = max(abs(x - xlo), abs(x - xhi)), max(abs(y - ylo), abs(y - yhi))
            assert dx * dx + dy * dy < 1
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("certificate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--input", type=Path, default=Path(__file__).with_name("outer_domain_certificate.json"))
    args = ap.parse_args()
    started = time.monotonic()
    raw = args.certificate.read_bytes()
    r = json.loads(raw)
    src = args.input
    source_raw = src.read_bytes()
    assert hashlib.sha256(source_raw).hexdigest() == r["source_global_box"]["sha256"]
    inp = json.loads(source_raw)
    den = int(r["coordinate_denominator"])
    assert den == int(inp["disk_polygon_denominator"])
    assert r["tau"] == str(TAU)
    disk = [(int(x), int(y), "D") for x, y in inp["disk_polygon_points"]]
    assert all(x * x + y * y <= den * den // 4 for x, y, _ in disk)
    q = inner_triangle()
    qs = [(int(x * den), int(y * den), f"Q{i}") for i, (x, y) in enumerate(q)]
    q_by_tag = {tag: (x, y) for x, y, tag in qs}
    assert all(F(x, den) == q[i][0] and F(y, den) == q[i][1] for i, (x, y, _) in enumerate(qs))
    lo, hi = -21 * den // 100, 21 * den // 100
    threshold2 = (2 * TAU * den * den).numerator
    max_level = int(r["max_level"])
    n = 2 ** max_level
    fine = (hi - lo) // n
    assert fine * n == hi - lo
    covered = bytearray(n * n)
    pruned_checked = 0
    retained = 0
    pruned_area = 0
    retained_area = 0
    for leaf in r["leaves"]:
        lev = int(leaf["level"])
        x0, x1 = map(int, leaf["x_scaled"])
        y0, y1 = map(int, leaf["y_scaled"])
        assert 0 <= lev <= max_level and lo <= x0 < x1 <= hi and lo <= y0 < y1 <= hi
        assert (x0 - lo) % fine == (x1 - lo) % fine == (y0 - lo) % fine == (y1 - lo) % fine == 0
        ix0, ix1, iy0, iy1 = (x0 - lo) // fine, (x1 - lo) // fine, (y0 - lo) // fine, (y1 - lo) // fine
        for ix in range(ix0, ix1):
            for iy in range(iy0, iy1):
                k = ix * n + iy
                assert not covered[k], "terminal cells overlap"
                covered[k] = 1
        assert leaf["status"] in ("prune", "retain")
        if leaf["status"] != "prune":
            retained_area += (x1-x0)*(y1-y0)
            retained += 1
            continue
        assert leaf["reason"] == "strict_inner_fan_area_exceeds_tau"
        assert (x0 + x1) % 2 == (y0 + y1) % 2 == 0
        cx, cy = (x0 + x1) // 2, (y0 + y1) // 2
        h = hull(disk + [(x + cx, y + cy, z) for x, y, z in qs])
        assert all(((q_by_tag[tag][0] + cx, q_by_tag[tag][1] + cy) if tag != "D" else (x, y)) == (x, y)
                   for x, y, tag in h)
        assert len(h) == leaf["fan_vertices_at_cell_center"]
        values = []
        for tx, ty in ((x0, y0), (x0, y1), (x1, y0), (x1, y1)):
            # h was formed at the centre.  Recover Q from its tag, instead
            # of adding a corner translation to an already-centred vertex.
            v = [(q_by_tag[tag][0] + tx, q_by_tag[tag][1] + ty) if tag != "D" else (x, y)
                 for x, y, tag in h]
            assert all(det(a, b) > 0 for a, b in zip(v, v[1:] + v[:1]))
            values.append(sum(det(a, b) for a, b in zip(v, v[1:] + v[:1])))
        assert values[0] + values[3] == values[1] + values[2]
        assert min(values) == int(leaf["minimum_twice_fan_area_scaled"])
        assert min(values) > threshold2
        pruned_area += (x1-x0)*(y1-y0)
        pruned_checked += 1
    assert all(covered), "terminal cells do not cover the full B3 translation square"
    assert pruned_area + retained_area == (hi-lo)**2
    ratios = r["translation_plane_area"]["fraction_of_original_translation_square"]
    assert F(ratios["pruned"]) == F(pruned_area,(hi-lo)**2)
    assert F(ratios["retained"]) == F(retained_area,(hi-lo)**2)
    assert r["terminal_leaf_counts"] == {"pruned":pruned_checked,"retained":retained,"all":len(r["leaves"])}
    result = {
        "status": "verified",
        "pruned_fraction_recomputed": str(F(pruned_area,(hi-lo)**2)),
        "retained_fraction_recomputed": str(F(retained_area,(hi-lo)**2)),
        "certificate": str(args.certificate),
        "certificate_sha256": hashlib.sha256(raw).hexdigest(),
        "checked_pruned_leaves": pruned_checked,
        "retained_leaves_preserved": retained,
        "partition": f"exact disjoint cover of 2^{max_level} by 2^{max_level} fine cells",
        "scope": "Verifies only the reported B3-plane prunings; retained leaves and other F8 parameters remain unresolved.",
        "python": sys.version,
        "elapsed_seconds": time.monotonic() - started,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
