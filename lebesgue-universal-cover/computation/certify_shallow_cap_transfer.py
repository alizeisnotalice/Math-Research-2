"""Certify the shallow cap comparison using the global five-body minimum."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
from zipfile import ZipFile

from certify_polygon_083731 import HERE, LAMBDA, ar, arb, ctx, diameter_squared
from tangent_polygons_finite_cover import generate


def read(name):
    return json.loads((HERE / name).read_bytes())


def certify(output):
    ctx.prec = 160
    seed = read("tomar-083740-seed.json")
    goal = F(seed["target_integer"], 2**57)
    assert F(".83740") < goal < F(".83740000000000001")
    assert all(F(row["integer_radius"], 2**44) < F(".069") and
               F(row["integer_area_lower"], 2**57) > goal for row in seed["radii"])
    roots = read("tomar-083740-retained.json")
    keys = {(int(k), v) for k, entries in roots.items() for v in entries}
    ordinary = read("finite-point-lift-083740.json")["rows"]
    exceptional = {(1205, p) for p in (139716, 151617, 139730, 152529)}
    assert len(keys) == 1940 and len(ordinary) == 1936
    assert {tuple(row[:2]) for row in ordinary} == keys - exceptional
    four = read("tomar-four-083740-margin-macbook.json")
    assert four["complete_four"] and F(four["target"]) == F(".837403")
    with ZipFile(HERE / "tomar-four-083740-margin.zip") as archive, ZipFile(HERE / "five-selected-083740-proofkit.zip") as previous:
        for part, pos in exceptional:
            stem = f"expected/part{part}-pos{pos}"
            raw_record = archive.read(stem+".json")
            record = json.loads(raw_record)
            assert sha256(raw_record).hexdigest() == four["leaves"][str(pos)]["record_sha256"]
            original = json.loads(previous.read(stem+".json"))
            for field in ("upstream", "original_tree_sha256", "original_header_sha256", "root"):
                assert record[field] == original[field]
            assert (record["part"], record["pos"]) == (part, pos)
            assert F(record["threshold"]) == F(".837403") and not record["result"]["remaining"]
            assert sha256(archive.read(stem+".bin")).hexdigest() == four["leaves"][str(pos)]["trace_sha256"]
    five = read("five-global-minimum-single-worker-replay.json")
    assert five["complete"] and arb(five["minimum_interval"]) > ar(F(".837295"))
    s, r, R = LAMBDA, F(1, 2), F(".703")
    eps, rho = 1-s, F(".000000735")
    circ = 1/(2*(arb.pi()/6).cos())
    assert circ < ar(F(".6"))
    displacement = F(".6")*eps
    steiner = F(22, 7)*(2*R*displacement+displacement**2)
    assert F(".837403")-steiner > goal
    radius = ar(s*r)
    radius_area = (arb.pi()-(radius/ar(R)).acos())*radius**2 + radius*(ar(R)**2-radius**2).sqrt()
    assert radius_area > ar(goal)
    height = F(".013")
    excess_area = goal-s**2*F(".837295")
    assert s*r*height**2/R > excess_area
    source = read("finite-point-leaf-input.json")
    points = source["points"]
    denominator = source["denominator"]
    qradius = F(".528619")
    assert all((F(x,denominator)+F(".01007"))**2+F(y,denominator)**2 < qradius**2 for x,y in points)
    m, a = rho*qradius/eps, r-F(".069")
    assert m < a < r and R**2 < 2*(s*r)**2
    v = (1+s)*r*r+m*m*eps
    derivative = 1+s-(2*m*r+v)/(4*r*r)
    right = (1+s)*height+max(0,(v-2*m*(r+height))/(2*r+height))
    margin = 4*(a-m)*s*s-right
    assert derivative > 0 and margin > F(".0017")

    diameter = diameter_squared(source)
    grid = 10**18
    scaled_point = F(isqrt(grid**2*diameter.denominator//diameter.numerator), grid)
    factor = scaled_point/(1-rho)
    groups, budget = generate((s*factor)**2, F(1))
    assert F(budget["common_dilation"])*s*factor < 1
    groups[-1] = dict(body="P30", denominator=denominator*scaled_point.denominator,
                      points=[[scaled_point.numerator*x,scaled_point.numerator*y] for x,y in points])
    assert diameter_squared(groups[-1]) <= 1
    floor = factor**2*goal
    assert floor > F(".8374012322")
    raw = (json.dumps(groups,separators=(",",":"))+"\n").encode()
    output.write_bytes(raw)
    return dict(goal=str(goal), rho=str(rho), cap_height=str(height), five_area_excess=str(excess_area),
                four_loss=str(steiner), comparison_margin=str(margin), radius_area=str(radius_area),
                common_scale=str(factor), area_floor=str(floor), counts=budget["counts"],
                squared_diameter_numerators=budget["squared_diameter_numerators"],
                polygon_denominator=budget["denominator"], P30_squared_diameter=str(diameter_squared(groups[-1])),
                sha256=sha256(raw).hexdigest())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    print(json.dumps(certify(parser.parse_args().output),indent=2))
