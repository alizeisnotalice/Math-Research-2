"""Construct a rational certificate at a target within the R152 contracted-family bound."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from certify_polygon_083731 import LAMBDA, diameter_squared
from tangent_polygons_finite_cover import HERE, ctx, generate


def transfer(target):
    binding = (HERE/"joint-contraction-audit.json").read_bytes()
    audit = json.loads(binding)
    assert F(audit["rho"]) == F(".000000555")
    assert sha256((HERE/"joint-contraction-points.json").read_bytes()).hexdigest() == audit["point_sha256"]
    goal = F(audit["threshold"])
    fixed = json.loads((HERE/"finite-point-leaf-input.json").read_bytes())
    original = {"points":fixed["points"], "denominator":fixed["denominator"]}
    ratio = target/goal
    assert 0 < ratio
    grid = 10**18
    factor = F(isqrt(grid**2*ratio.numerator//ratio.denominator)+1,grid)
    sixth_scale = factor*(1-F(audit["rho"]))
    assert sixth_scale**2*diameter_squared(original) <= 1
    width = factor*LAMBDA
    assert 0 < width < 1
    groups, budget = generate(width**2,F(1))
    assert F(budget["common_dilation"])*width < 1
    groups[-1] = {"body":"P30","points":[[x*sixth_scale.numerator,y*sixth_scale.numerator]
                     for x,y in original["points"]],"denominator":original["denominator"]*sixth_scale.denominator}
    diameters = [diameter_squared(g) for g in groups]
    assert all(d <= 1 for d in diameters)
    floor = goal*factor**2
    assert floor > target
    return groups,dict(target=str(target), area_floor=str(floor), factor=str(factor),
                       curved_core_scale=str(width), counts=budget["counts"], points=budget["total"],
                       squared_diameters=list(map(str,diameters)), core_transfer_sha256=sha256(binding).hexdigest())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target",type=F)
    parser.add_argument("coordinates",type=Path)
    args = parser.parse_args()
    ctx.prec = 160
    groups, audit = transfer(args.target)
    raw = (json.dumps(groups,separators=(",",":"))+"\n").encode()
    args.coordinates.write_bytes(raw)
    audit["coordinate_sha256"] = sha256(raw).hexdigest()
    print(json.dumps(audit,indent=2))
