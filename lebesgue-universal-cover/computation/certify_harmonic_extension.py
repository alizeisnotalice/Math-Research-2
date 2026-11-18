"""Replay all-angle lower bounds, a pointwise failure, and one free eight-body upper bound."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from math import lcm
from pathlib import Path
import platform

import flint
from flint import arb, ctx

import certify_irregular_hex
import certify_nine
from certify_irregular_hex import NS, RHO_PAD, determinant, support
from certify_nine import FOURIER_TERMS, PI_UPPER, S, append_point, ar, arc_body, arithmetic, cross, fourier_body, hull

HERE = Path(__file__).resolve().parent
TERMS = ((3, F(1, 80)), (5, F(13, 800)))
FAR_LIMIT = F(1, 4)
INNER_RADIUS = F(47, 100)
ROTATION_LIPSCHITZ = F(19, 160)
LOWER_TARGET = F("0.8386")
POINT_BUDGET = F("0.83809903")
POINT_UPPER_TARGET = F("0.83799238")
FREE_EIGHT_TARGET = F("0.837723081")


def disk_points(points, count, pi):
    for j in range(count):
        angle = 2*pi*j/count
        append_point(points, angle.cos()/2, angle.sin()/2)


def host_polygon(spec, k, pi):
    """Enclose samples at the six Round 65 root intervals, for any choice of roots."""
    d1, d2 = map(ar, spec["angle_perturbations"])
    angles = (arb(0), 2*pi/3 + d1, 4*pi/3 + d2)
    normals = [(a.cos(), a.sin()) for a in angles]
    weights = [determinant(normals[(j+1)%3], normals[(j+2)%3]) for j in range(3)]
    assert all(w > 0 for w in weights) and normals[1][1] > 0
    assert len(spec["rotation_root_proposals"]) == 6

    points = []
    disk_points(points, 8*k, pi)
    for j, center in enumerate(map(F, spec["rotation_root_proposals"])):
        def equation(phi):
            return sum(w*(support(j, a-phi, pi)-ar(F(1, 2)))
                       for w, a in zip(weights, angles, strict=True))

        lo, hi = ar(center - RHO_PAD), ar(center + RHO_PAD)
        left, right = equation(lo), equation(hi)
        assert (left < 0 < right) or (right < 0 < left)
        rho = arb.union(lo, hi)
        tx = ar(F(1, 2)) - support(j, -rho, pi)
        ty = (ar(F(1, 2)) - support(j, angles[1]-rho, pi) - tx*normals[1][0])/normals[1][1]
        if j < 4:
            arc_body(points, k, pi, NS[j], rho, tx, ty)
        else:
            fourier_body(points, (8, 16)[j-4]*k, pi, FOURIER_TERMS[j-4], rho, tx, ty)
    return hull(points)


def marked_hull(points):
    """Andrew hull retaining the fixed or translated label on each selected point."""
    unique = []
    for point in sorted(points):
        if not unique or point[:2] != unique[-1][:2]:
            unique.append(point)

    def half(sequence):
        out = []
        for point in sequence:
            while len(out) > 1 and cross(out[-2], out[-1], point) <= 0:
                out.pop()
            out.append(point)
        return out[:-1]

    poly = half(unique) + half(reversed(unique))
    assert len(poly) > 2
    assert all(cross(a, b, c) > 0 for a, b, c in
               zip(poly, poly[1:]+poly[:1], poly[2:]+poly[:2], strict=True))
    return poly


def angular_lower(old, material, left, right, witnesses, turn_den, shift_den, pi):
    """Recompute supporting affine cuts and a weighted lower bound for one angle cell."""
    a, b = F(left, turn_den), F(right, turn_den)
    angle = pi*ar(a+b)
    co, si = angle.cos(), angle.sin()
    rotated = []
    for px, py in material:
        x, y = ar(F(px, S)), ar(F(py, S))
        append_point(rotated, x*co-y*si, x*si+y*co)

    error = F(4, S) + ROTATION_LIPSCHITZ*PI_UPPER*(b-a)
    shrink = INNER_RADIUS/(INNER_RADIUS+error)
    old_den, new_den = S+4, S*shrink.denominator
    scale = lcm(old_den, new_den, shift_den)
    fixed = [(x*(scale//old_den), y*(scale//old_den), 0) for x, y in old]
    factor = shrink.numerator*(scale//new_den)
    moved = [(x*factor, y*factor) for x, y in hull(rotated)]

    value = gx = gy = weight_sum = 0
    for tx, ty, weight in witnesses:
        assert all(type(v) is int for v in (tx, ty, weight))
        assert weight > 0
        x, y = tx*(scale//shift_den), ty*(scale//shift_den)
        poly = marked_hull(fixed+[(u+x, v+y, 1) for u, v in moved])
        twice = dx = dy = 0
        for (px, py, label), (qx, qy, next_label) in zip(poly, poly[1:]+poly[:1], strict=True):
            twice += px*qy-py*qx
            dx += label*qy-next_label*py
            dy += next_label*px-label*qx
        value += weight*(twice-dx*x-dy*y)
        gx += weight*dx
        gy += weight*dy
        weight_sum += weight
    assert weight_sum > 0
    return F(value, 2*scale**2*weight_sum) - FAR_LIMIT*F(abs(gx)+abs(gy), 2*scale*weight_sum)


def full_lower(specs, inputs, pi):
    s = ar(FAR_LIMIT + INNER_RADIUS)
    far_area = (pi-(1/(2*s)).acos())/4+(s**2-arb(1)/4).sqrt()/2
    assert far_area > ar(LOWER_TARGET)

    material = []
    fourier_body(material, 384, pi, TERMS, arb(0), arb(0), arb(0))
    turn_den = inputs["turn_denominator"]
    shift_den = inputs["translation_denominator"]
    assert type(turn_den) is int and turn_den > 0
    assert type(shift_den) is int and shift_den > 0
    rows = []
    for spec, cells in zip(specs, inputs["angular_witnesses"], strict=True):
        old = host_polygon(spec, 96, pi)
        end = 0
        bounds = []
        for stop, witnesses in cells:
            assert type(stop) is int and end < stop <= turn_den
            bound = angular_lower(old, material, end, stop, witnesses, turn_den, shift_den, pi)
            assert bound > LOWER_TARGET
            bounds.append(bound)
            end = stop
        assert end == turn_den
        rows.append({
            "angle_perturbations": spec["angle_perturbations"],
            "angular_cells": len(bounds),
            "host_inner_polygon_vertices": len(old),
            "lower": str(min(bounds)),
            "display": float(min(bounds)),
        })
    return rows


def single_point_upper(spec, pose, pi):
    """Verify one entire extra copy in the pointwise area-budget region."""
    k, count = 512, 2048
    old = host_polygon(spec, k, pi)
    old_error = F(5, 256)*PI_UPPER**2/k**2 + F(2, S)
    new = []
    fourier_body(new, count, pi, TERMS, *map(ar, pose))
    add_error = F(3, 8)*(2*PI_UPPER/count)**2 + F(2, S)
    scale = lcm(S, add_error.denominator)
    outer_offset = int(add_error*scale)

    fixed = [(x*(scale//S), y*(scale//S)) for x, y in old]
    outer = hull([(x*(scale//S)+a*outer_offset, y*(scale//S)+b*outer_offset)
                  for x, y in hull(new) for a in (-1, 1) for b in (-1, 1)])
    assert all(F(x*x+y*y, scale*scale) < F(9, 25) for x, y in fixed+outer)
    assert 2*PI_UPPER*F(3, 5) < 4

    twice, edges = 0, []
    for (x, y), (u, v) in zip(fixed, fixed[1:]+fixed[:1], strict=True):
        a, b = v-y, x-u
        edges.append((a, b, a*x+b*y))
        twice += x*v-y*u
    increment = max(sum(max(a*x+b*y-c, 0) for a, b, c in edges) for x, y in outer)
    upper = F(twice+increment, 2*scale**2) + 4*old_error + PI_UPPER*old_error**2
    assert upper < POINT_UPPER_TARGET < POINT_BUDGET
    return {"outer_vertices": len(outer), "upper": str(upper), "display": float(upper)}


def free_eight_upper(pose, pi):
    k = 2048
    assert len(pose) == 20
    t = list(map(ar, pose))
    points = []
    disk_points(points, 8*k, pi)
    arc_body(points, k, pi, 3, arb(0), *t[:2])
    for n, j in ((5, 2), (7, 5), (9, 8)):
        arc_body(points, k, pi, n, *t[j:j+3])
    for (coeffs, j), segments in zip(((TERMS, 11), (FOURIER_TERMS[0], 14), (FOURIER_TERMS[1], 17)),
                                      (16, 8, 16), strict=True):
        fourier_body(points, segments*k, pi, coeffs, *t[j:j+3])
    poly = hull(points)
    result = arithmetic(poly, k)
    assert F(result["upper"]) < FREE_EIGHT_TARGET < POINT_BUDGET
    return {"sample_count": len(points), "polygon_vertices": len(poly), **result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="also compare the archived exact result")
    args = parser.parse_args()
    ctx.prec = 160
    pi = arb.pi()
    assert pi < ar(PI_UPPER)
    amplitude = sum((k*k-1)*abs(a) for k, a in TERMS)
    derivative = sum(k*(k*k-1)*abs(a) for k, a in TERMS)
    assert F(1, 2)-amplitude == F(1, 100)
    assert (F(1, 2)+amplitude)**2+derivative**2 < 9
    assert F(3, 8*8**2) < F(5, 256)
    assert F(1, 2)-sum(abs(a) for _, a in TERMS) > INNER_RADIUS
    assert sum(k*abs(a) for k, a in TERMS) == ROTATION_LIPSCHITZ

    raw = (HERE/"irregular_inputs.json").read_bytes()
    raw_inputs = (HERE/"harmonic_extension_inputs.json").read_bytes()
    specs = json.loads(raw)
    inputs = json.loads(raw_inputs)
    assert sha256(raw).hexdigest() == inputs["round65_input_sha256"] and len(specs) == 3
    host = inputs["single_point_host"]
    assert type(host) is int and host == 2

    old_lipschitz = F(29, 100)
    assert old_lipschitz**2 > F(1, 12)
    assert all(sum(k*abs(a) for k, a in terms) < old_lipschitz for terms in FOURIER_TERMS)
    rho = (F(99, 70)+old_lipschitz)/12000
    assert F(99, 70)**2 > 2
    neighborhood = (LOWER_TARGET-PI_UPPER*rho**2)/(1+4*rho)
    assert neighborhood > F("0.83812")
    result = {
        "precision_bits": ctx.prec,
        "dependencies_sha256": {
            "irregular_inputs.json": sha256(raw).hexdigest(),
            "harmonic_extension_inputs.json": sha256(raw_inputs).hexdigest(),
            "certify_irregular_hex.py": sha256(Path(certify_irregular_hex.__file__).read_bytes()).hexdigest(),
            "certify_nine.py": sha256(Path(certify_nine.__file__).read_bytes()).hexdigest(),
            "certify_harmonic_extension.py": sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "exact_results": {
            "fixed_hosts_all_placements": full_lower(specs, inputs, pi),
            "fixed_host_common_strict_threshold": str(LOWER_TARGET),
            "host_neighborhood_halfwidth": "1/12000",
            "host_neighborhood_lower": str(neighborhood),
            "single_point_budget": str(POINT_BUDGET),
            "single_point": single_point_upper(specs[host], inputs["single_point_pose"], pi),
            "freely_rearranged_eight": free_eight_upper(inputs["free_eight_pose"], pi),
        },
        "runtime": {"python": platform.python_version(), "python_flint": flint.__version__},
    }
    if args.check:
        archived = json.loads((HERE/"harmonic_extension_results.json").read_bytes())
        assert {k: v for k, v in archived.items() if k != "runtime"} == {
            k: v for k, v in result.items() if k != "runtime"
        }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
