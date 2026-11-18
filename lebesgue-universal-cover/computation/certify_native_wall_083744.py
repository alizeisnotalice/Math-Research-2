"""Verify the local translation bound used by the positive P30 terminals at area .83744."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys


def load(path, expected=None):
    raw = path.read_bytes()
    if expected:
        assert sha256(raw).hexdigest() == expected
    return json.loads(raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("computation", "tree", "table32", "table64", "work"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    assert __debug__
    comp = args.computation.resolve()
    sys.path.insert(0, str(comp))
    from certify_polygon_083731 import ORDERS, ar, arb, ctx, lower_radius
    from certify_refined_material_contraction import verify as check_material
    from certify_shallow_cap_transfer import diameter_squared
    from certify_tomar_six_high import SourceTree
    from tangent_polygons_finite_cover import generate

    ctx.prec = 160
    args.work.mkdir(parents=True, exist_ok=False)
    tree = SourceTree(args.tree)
    terminal = load(comp / "finite-point-lift-083744.json", "af6e8b2bff8b52408351fa1e0558358311b4ac64ff6fb0c0f3dc619243a54a6e")
    previous = load(comp / "signed-shallow-cap-083744-replay.json")
    assert previous["complete"] and previous["ordinary"] == 2833 and previous["pure"] == 124
    assert sha256((comp / "finite-point-lift-083744.json").read_bytes()).hexdigest() == previous["terminal_sha256"]
    goal = F(previous["goal"])
    assert goal == F(120687823102244742, 2**57)
    dense = load(comp / terminal["dense"]["filename"], terminal["dense"]["sha256"])
    scale, radius, outer = F(".999999006"), F(1, 2), F(".703")
    assert F(dense["scale"]) == scale and len(dense["receipts"]) == 6
    material = check_material(comp, args.tree, args.table32, args.table64)
    assert material == load(comp / "refined-material-contraction.json")
    assert F(material["curve_scale"]) == scale
    ordinary = terminal["rows"]
    assert len({tuple(row[:2]) for row in ordinary}) == len(ordinary)
    keys = {tuple(row[:2]) for row in ordinary if row[5] > 0}
    assert keys.isdisjoint(map(tuple, dense["original_leaves"]))
    keys.update(map(tuple, dense["original_leaves"]))
    assert len(keys) == 1636
    max_square, achiever = F(0), None
    for part, pos in sorted(keys):
        box = tree.leaf(part, pos)
        square = sum(max(abs(lo), abs(hi))**2 for lo, hi in box[:2])
        if square > max_square:
            max_square, achiever = square, [part, pos]
    wall = F(".47184")
    assert max_square < (radius-wall)**2
    pi = arb.pi()
    assert all(lower_radius(n, pi) > ar(wall) for n in (5, 7, 9))
    sr = scale*radius
    assert pi*(ar(sr)**2) - ar(sr)**2*(ar(sr/outer)).acos() + ar(sr)*(ar(outer**2-sr**2)).sqrt() > ar(goal)
    assert outer**2 < 2*sr**2
    for offset, n in enumerate(ORDERS[1:]):
        k = 3*offset
        displacement = min(F(tree.root["hi"][k]), F(tree.root["hi"][k+1]),
                           -F(tree.root["lo"][k]), -F(tree.root["lo"][k+1]))
        inner = ar(scale)*lower_radius(n, pi)
        assert pi*(ar(sr)**2 + inner**2)/2 + ar(displacement)*(ar(sr)+inner) > ar(goal)
    seed = load(comp / "tomar-083744-seed.json")
    assert len(seed["radii"]) == 4 and all(F(row["integer_area_lower"], 2**57) > goal for row in seed["radii"])
    delta = F(".6")*(1-scale)
    loss = F(22, 7)*(2*outer*delta+delta**2)
    assert len(previous["high_rows"]) == 124 and all(F(row[2])-loss > goal for row in previous["high_rows"])

    rho, h, q, denominator = F(".000000860"), F(".0121"), F(".528619"), F(".495")
    five = load(comp / "five-global-minimum-single-worker-replay.json")
    assert five["complete"] and arb(five["minimum_interval"]) > ar(F(".837295"))
    assert denominator**2 > outer**2-sr**2
    allowed = scale**2*F(".837295")+sr*h**2/denominator
    assert allowed > goal
    cloud = load(comp / "finite-point-leaf-input.json", terminal["point_sha256"])
    n = cloud["denominator"]
    assert all((F(x,n)+F(".01007"))**2+F(y,n)**2 < q*q for x, y in cloud["points"])
    m = rho*q/(1-scale)
    v = (1+scale)*radius**2 + (1-scale)*m**2
    right = (1+scale)*h + max(0, (v-2*m*(radius+h))/(2*radius+h))
    margin = 4*(wall-m)*scale**2-right
    assert m < wall < radius
    assert 1+scale-(2*m*radius+v)/(4*radius*radius) > 0 and margin > F(".002536")

    # The previous fixed point coordinates bind both exact independent dilation factors.
    fixed = load(comp / "shallow-cap-points.json", previous["coordinate_sha256"])
    assert len(fixed) == 6 and all(g["denominator"] == 1000000229098 for g in fixed[:5])
    # The coordinate representation uses numerator/denominator of the normalized P30 scale.
    from math import isqrt
    diam = diameter_squared(cloud)
    t = F(isqrt(10**36*diam.denominator//diam.numerator), 10**18)
    assert fixed[5]["denominator"] == t.denominator*n
    assert fixed[5]["points"] == [[t.numerator*x, t.numerator*y] for x,y in cloud["points"]]
    factor = t/(1-rho)
    fixed_factor = min(factor, F(10**12-5, 1000000229098)/scale)
    assert sum(len(z["points"]) for z in fixed) == 13368
    assert fixed_factor**2*goal > F(".83744128111")
    groups, budget = generate((scale*factor)**2, F(1))
    assert F(budget["common_dilation"])*scale*factor < 1
    groups[-1] = fixed[-1]
    assert diameter_squared(groups[-1]) <= 1 and factor**2*goal > F(".83744144166")
    coordinates = (json.dumps(groups,separators=(",",":"))+"\n").encode()
    (args.work / "coordinates.json").write_bytes(coordinates)
    result = dict(complete=True, native_point_leaves=len(keys), max_b3_norm_squared=str(max_square),
                  attaining_leaf=achiever, minimum_wall=str(wall), positive_margin=str(margin),
                  radius_denominator=str(denominator), maximum_cap_goal=str(allowed),
                  old_five_loss=str(loss), scale=str(scale), point_shrink=str(rho), goal=str(goal),
                  factor=str(factor), fixed_factor=str(fixed_factor), fixed_lower=str(fixed_factor**2*goal),
                  area_lower=str(factor**2*goal), counts=[len(g["points"]) for g in groups],
                  fixed_coordinate_sha256=previous["coordinate_sha256"],
                  coordinate_sha256=sha256(coordinates).hexdigest(),
                  ordinary_input_sha256=previous["terminal_sha256"], dense_sha256=terminal["dense"]["sha256"])
    (args.work/"audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
