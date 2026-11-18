"""Check two physical Reuleaux hulls and quantify perturbations of the centered disk/triangle model."""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx


def main(output, subdivisions):
    comp = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(comp))
    from certify_nine import ar, append_point, arc_body, arithmetic, hull

    ctx.prec = 192
    pi = arb.pi()
    half = ar(Fraction(1, 2))
    root = arb(3).sqrt()
    radius = half / (pi/10).cos()
    excess = radius*(pi/30).cos() - half
    base = []
    for j in range(8*subdivisions):
        sine, cosine = (pi*j/(4*subdivisions)).sin_cos()
        append_point(base, half*cosine, half*sine)
    for j in range(3):
        sine, cosine = (pi*(1+4*j)/6).sin_cos()
        append_point(base, cosine/root, sine/root)
    poses = {"full_optimum": (pi/30, -excess, excess/root),
             "vertex_optimum": (pi/10, arb(0), 2*excess/root)}
    cases = {}
    for label, (angle, tx, ty) in poses.items():
        points = base.copy()
        arc_body(points, subdivisions, pi, 5, angle, tx, ty)
        poly = hull(points)
        result = arithmetic(poly, subdivisions)
        cases[label] = {key: result[key] for key in ("lower", "upper", "epsilon")}
        cases[label]["hull_vertices"] = len(poly)
    exact = json.loads(Path(__file__).with_name("area.json").read_text())["full_minimum"]
    low, high = map(Fraction, exact)
    opt = cases["full_optimum"]
    at_vertex_pose = cases["vertex_optimum"]
    assert Fraction(opt["lower"]) < low < high < Fraction(opt["upper"]) < Fraction(at_vertex_pose["lower"])

    # Each of the five boundary arcs subtends pi/5 at its center and has unit radius.
    # With a disk/triangle center mismatch bounded by d, translate the two hosts
    # by half the offset relative to each other before taking Hausdorff distance.
    targets = []
    for threshold, center_gap in [("0.834", "0"), ("0.835", "0"), ("0.835", "0.00007"),
                                  ("0.83510", "0.00001"), ("0.83513", "0")]:
        target, gap = ar(threshold), ar(center_gap)
        for count in range(1, 1000):
            sagitta = 1 - (pi/(10*count)).cos()
            error = sagitta.max(gap/2)
            floor = ar(low)/(1 + 2*error)**2
            if floor > target:
                assert count == 1 or ar(high)/(1 + 2*(1-(pi/(10*(count-1))).cos()).max(gap/2))**2 < target
                targets.append({"strict_target": threshold, "center_gap": center_gap,
                                "segments_per_arc": count, "input_points": 5*count,
                                "certified_floor": str(floor)})
                break
        else:
            raise ValueError(threshold)
    output.write_text(json.dumps({"subdivisions": subdivisions, "cases": cases, "finite_targets": targets},
                                 indent=2) + "\n")
    print([(k, float(Fraction(v["lower"])), float(Fraction(v["upper"]))) for k,v in cases.items()])
    print([(v["strict_target"], v["center_gap"], v["input_points"]) for v in targets])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--subdivisions",type=int,default=8192)
    args = parser.parse_args()
    main(args.output,args.subdivisions)
