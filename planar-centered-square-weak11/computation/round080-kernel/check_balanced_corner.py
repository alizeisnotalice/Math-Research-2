#!/usr/bin/env python3
"""R80 exact raw-source cap qualification and scalar analytic-area certificate.

Does not compute a rectangle union numerically or call an area producer.
"""
import argparse
from fractions import Fraction as Q
import hashlib
from itertools import combinations_with_replacement, product
import json
from math import sqrt
import os
from pathlib import Path
import platform
import sys
import time


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def plain(v):
    if isinstance(v, Q):
        return str(v)
    if isinstance(v, dict):
        return {str(k): plain(x) for k, x in v.items()}
    if isinstance(v, (tuple, list)):
        return [plain(x) for x in v]
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    started = time.perf_counter()
    data = json.loads(args.input.read_text())
    pts = [tuple(map(Q, x)) for x in data["positions"]]
    masses = list(map(Q, data["masses"]))
    assert data["threshold"] == "1" and data["S_indices"] == [0, 1, 2, 3]
    w, e, d, D, v = Q(999, 1000), Q(1, 1000), Q(71, 200), Q(171, 200), Q(1, 4)
    assert pts == [(Q(0), Q(0)), (w, Q(0)), (Q(0), w), (w, w), (-d, Q(0)), (Q(0), -D)]
    assert masses == [v]*5+[Q(499, 1000)]
    assert list(map(Q, data["corner_crop"])) == [-Q(1, 2), -e, -Q(1, 2), -e]
    N = len(pts)

    def ids(mask):
        return [i for i in range(N) if mask & (1 << i)]

    def bbox(mask):
        I = ids(mask)
        return [min(pts[i][0] for i in I), max(pts[i][0] for i in I),
                min(pts[i][1] for i in I), max(pts[i][1] for i in I)]

    def closure(box):
        a, b, c, f = box
        return sum(1 << i for i, (x, y) in enumerate(pts) if a <= x <= b and c <= y <= f)

    xx, yy = sorted(set(x for x, y in pts)), sorted(set(y for x, y in pts))
    box_masks = [closure((a, b, c, f))
                 for (a, b), (c, f) in product(combinations_with_replacement(xx, 2),
                                               combinations_with_replacement(yy, 2))]
    assert len(box_masks) == 36
    closures = set(box_masks)-{0}
    subset_closures = {closure(bbox(mask)) for mask in range(1, 1 << N)}
    assert closures == subset_closures and len(closures) == 23
    rows = []
    for mask in sorted(closures):
        a, b, c, f = box = bbox(mask)
        M = sum((masses[i] for i in ids(mask)), Q(0))
        span = max(b-a, f-c)
        slack = max(Q(1), span*span)-M
        assert slack >= 0
        rows.append({"mask": mask, "bbox": box, "mass": M, "span_squared": span*span,
                     "cap_slack": slack, "eligible": M > span*span})
    eligible = [row for row in rows if row["eligible"]]
    maximum = max(row["mass"] for row in eligible)
    maximum_masks = [row["mask"] for row in eligible if row["mass"] == maximum]
    assert maximum == 1 and maximum_masks == [15]
    bymask = {row["mask"]: row for row in rows}
    assert all(bymask[mask]["eligible"] for mask in [15, 17, 33])
    assert 15 | 17 | 33 == (1 << N)-1
    # The three eligible edges meet at index 0, hence one component.
    assert all(mask & 1 for mask in [15, 17, 33])
    total = sum(masses)
    a, b, c, f = bbox((1 << N)-1)
    span = max(b-a, f-c)
    nonterminal_gap = span*span-Q(35, 33)*total
    assert total == Q(1749, 1000) and span == Q(927, 500)
    assert nonterminal_gap == Q(395579, 250000) > 0

    bounds = {key: {k: Q(value) for k, value in row.items()}
              for key, row in data["radical_bounds"].items()}
    for row in bounds.values():
        assert 0 < row["lower"] < row["upper"]
        assert row["lower"]**2 < row["squared"] < row["upper"]**2
    al, au = bounds["a"]["lower"], bounds["a"]["upper"]
    Ul, Uu = bounds["U"]["lower"], bounds["U"]["upper"]
    Vl, Vu = bounds["V"]["lower"], bounds["V"]["upper"]
    assert e < d-au and d-al < v < al and au < d < Vl < Vu < Ul < Uu < Q(1, 2)
    assert Q(499, 4000) < (D-Q(1, 2))**2  # residual singleton 2 misses crop
    assert D-Uu > au and D-Vu > au and D < 2*Vl
    assert d < Vl  # both A2 and pair residual have X crop [e,V]
    assert Q(1, 4)+Q(499, 1000) == Q(3, 4)-Q(1, 1000) == Q(749, 1000)

    # Exact symbolic polynomials in a,U,V, reduced using their rational squares.
    squares = [bounds[key]["squared"] for key in ["a", "U", "V"]]

    def const(q):
        return {(0, 0, 0): Q(q)}

    def add(*ps):
        out = {}
        for p in ps:
            for power, coefficient in p.items():
                out[power] = out.get(power, Q(0))+coefficient
        return {power: q for power, q in out.items() if q}

    def mul(p, q):
        out = {}
        for power1, coef1 in p.items():
            for power2, coef2 in q.items():
                power = tuple(x+y for x, y in zip(power1, power2))
                coef = coef1*coef2
                for i, exponent in enumerate(power):
                    coef *= squares[i]**(exponent//2)
                reduced = tuple(exponent % 2 for exponent in power)
                out[reduced] = out.get(reduced, Q(0))+coef
        return {power: coefficient for power, coefficient in out.items() if coefficient}

    def scale(q, p):
        return mul(const(q), p)

    A, U, V = [{tuple(int(i == j) for i in range(3)): Q(1)} for j in range(3)]
    minus = lambda p: scale(-1, p)
    # First cropped component: A1 minus R1 plus actual pure-old thin strip.
    first = add(mul(add(scale(2, A), const(-d)), add(A, const(-e))),
                scale(-(v-e), add(A, const(-d+v))),
                scale(v-e, add(const(d-e), minus(A))))
    first_expected = add(const(Q(182323, 500000)), scale(-Q(171, 200), A))
    assert first == first_expected
    joint = add(mul(add(U, const(-e)), add(scale(2, U), const(-D))),
                minus(mul(add(V, const(-e)), add(scale(2, V), const(-D)))))
    joint_expected = add(const(Q(1, 8)), scale(-Q(857, 1000), add(U, minus(V))))
    assert joint == joint_expected
    complete = add(first, joint)
    expected = add(const(Q(244823, 500000)), scale(-Q(171, 200), A),
                   scale(-Q(857, 1000), U), scale(Q(857, 1000), V))
    assert complete == expected
    first_lower = Q(182323, 500000)-Q(171, 200)*au
    joint_lower = Q(1, 8)-Q(857, 1000)*(Uu-Vl)
    lower = first_lower+joint_lower
    assert first_lower == Q(6235733, 100000000)
    assert joint_lower == Q(67557861, 1000000000)
    assert lower == Q(129915191, 1000000000) > Q(1, 8)
    gap = lower-Q(1, 8)
    assert gap == Q(4915191, 1000000000)
    approximate = float(Q(244823, 500000))-float(Q(171, 200))*sqrt(float(squares[0])) \
                  -float(Q(857, 1000))*(sqrt(float(squares[1]))-sqrt(float(squares[2])))
    result = {"status": "PASS", "scope": data["scope"],
              "global_qualification": {"coordinate_box_count": len(box_masks),
                "nonempty_subset_count": (1 << N)-1, "bbox_closure_count": len(closures),
                "closed_caps_all_pass": True, "bbox_closures": rows,
                "eligible_closure_count": len(eligible), "maximum_eligible_mass": maximum,
                "maximum_masks": maximum_masks, "unique_old_S": True,
                "connected": True, "connected_eligible_hyperedges": [15, 17, 33],
                "total_mass": total, "source_span": span, "nonterminal_gap": nonterminal_gap},
              "local_exact_geometry": {"full_local_source_count": 3, "full_label_count": 7,
                "complete_residual_label_count": 3, "old_local_mass": Q(1, 4),
                "radical_bounds": bounds, "A2_equals_pair_residual_on_crop": True,
                "old_only_thin_strip_included": True, "symbolic_area_coefficients_a_U_V": complete,
                "first_piece_strict_lower": first_lower, "joint_piece_strict_lower": joint_lower,
                "actual_SW_area_strict_lower": lower, "half_mass_comparison": Q(1, 8),
                "strict_surplus_lower": gap, "approximate_area_diagnostic_only": approximate},
              "receipt": {"pid": os.getpid(), "python": sys.version,
                "platform": platform.platform(), "actual_command": [sys.executable]+sys.argv,
                "working_directory": str(Path.cwd()), "input_sha256": sha(args.input),
                "script_sha256": sha(__file__), "wall_seconds": time.perf_counter()-started,
                "random_seed": None, "tolerance": "exact rational, no tolerance",
                "no_area_union_engine": True}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plain(result), indent=2)+"\n")
    print(json.dumps({"status": "PASS", "closures": len(closures),
                      "unique_maximum_mask": maximum_masks,
                      "actual_SW_strict_lower": str(lower), "strict_surplus_lower": str(gap)}))


if __name__ == "__main__":
    main()
