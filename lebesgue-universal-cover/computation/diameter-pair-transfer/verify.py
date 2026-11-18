"""Shorten separated diameter pairs, certify the transfer and construct six rational witnesses."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
import json
from math import isqrt
from pathlib import Path


def module(name, path):
    spec = spec_from_file_location(name, path)
    result = module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def upper_sqrt(q, digits=24):
    scale = 10**digits
    n = isqrt(q.numerator*scale**2 // q.denominator)
    if n*n*q.denominator < q.numerator*scale**2:
        n += 1
    return F(n, scale)


def run(comp, work):
    assert __debug__
    project = comp / "projected-cap-transfer"
    reference = json.loads((project / "verified.json").read_bytes())
    cap = module("projected_cap", project / "verify.py")
    assert cap.certify(comp, F(reference["point_shrink"])) == reference
    s, goal = map(F, (reference["scale"], reference["threshold"]))
    budget = (1-s)*F(reference["projection_rate"])
    points_file = comp / "finite-point-leaf-input.json"
    assert sha256(points_file.read_bytes()).hexdigest() == reference["sources"][points_file.name]
    source = json.loads(points_file.read_bytes())
    points, den = source["points"], source["denominator"]
    lengths = {(i, j): F((x-u)**2+(y-v)**2, den**2)
               for (i, (x,y)), (j, (u,v)) in combinations(enumerate(points), 2)}
    diameter = max(lengths.values())
    assert diameter == F(reference["point_diameter_squared"])
    match = [pair for pair, length in lengths.items() if length > F(".9999991")**2]
    partners = dict(match)
    partners.update((j, i) for i, j in match)
    assert len(match) == 15 and len(partners) == len(points) == 30
    secondary = max(length for pair, length in lengths.items() if pair not in match)
    du, nonmatch = upper_sqrt(diameter), upper_sqrt(secondary)
    assert nonmatch + 4*budget < F(".9999999") and F(".9999999")**2 < diameter
    rho = 2*budget/du
    # p' = (1-rho/2)*p + (rho/2)*p_partner with one integer denominator.
    u, v = (rho/2).as_integer_ratio()
    shifted = [[(v-u)*x+u*y for x,y in zip(point, points[partners[i]])]
               for i, point in enumerate(points)]
    maker = module("joint_generator", project / "generate_joint_cloud.py")
    geometry, _ = maker.load_geometry()
    moved_diameter, axis = geometry.farthest(shifted)
    assert F(moved_diameter, (v*den)**2) == diameter*(1-rho)**2
    for point, new in zip(points, shifted):
        assert sum((x*v-y)**2 for x,y in zip(point, new)) <= (v*den)**2*budget**2
    f2 = 1/(diameter*(1-rho)**2)
    assert s*s*f2 < 1
    work.mkdir(parents=True, exist_ok=False)
    maker.main(comp, work/"curve-data", rho, f2)
    base = json.loads((work/"curve-data/joint-receipt.json").read_bytes())
    groups = json.loads((work/"curve-data/joint-clouds.json").read_bytes())[:5]
    assert all(F(rec["inner_scale_squared"]) >= s*s*f2 for rec in base["group_certificates"][:5])
    groups.append(geometry.exact_unit_cloud("paired-P30", shifted, moved_diameter, axis))
    total = 0
    for group in groups:
        coords = group["points"]
        sq, _ = geometry.farthest(coords)
        assert sq == group["denominator"]**2
        total += len(coords)*(len(coords)-1)//2
    packed = geometry.packed(groups)
    (work/"points.json").write_bytes(packed)
    lower = goal*f2
    rounded = (lower*10**24).__floor__()
    audit = {
        "threshold": str(goal), "curve_scale": str(s), "source_point_sha256": sha256(points_file.read_bytes()).hexdigest(),
        "old_projection_budget": str(budget), "same_29_cap_slabs": reference["projection_slabs"],
        "matched_point_indices": match, "diameter_squared": str(diameter),
        "diameter_upper": str(du), "nonmatching_distance_upper": str(nonmatch),
        "common_contraction": str(rho), "modified_diameter_squared": str(diameter*(1-rho)**2),
        "factor_squared": str(f2), "lower": str(lower),
        "strict_decimal_lower": f"{rounded//10**24}.{rounded%10**24:024d}",
        "curve_containment": base["group_certificates"][:5],
        "counts": [len(g["points"]) for g in groups], "all_diameters_squared": [1]*6,
        "integer_point_pairs": total, "points_sha256": sha256(packed).hexdigest(),
    }
    assert F(audit["strict_decimal_lower"]) < lower
    (work/"verified.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(audit["strict_decimal_lower"], total, audit["counts"])


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("comp", type=Path)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    run(a.comp.resolve(), a.work.resolve())
