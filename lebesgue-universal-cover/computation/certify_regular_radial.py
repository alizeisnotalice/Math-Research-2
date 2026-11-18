"""Certify a uniform radial exclusion on translations of a regular Reuleaux polygon."""

import argparse
import json
import platform
from fractions import Fraction as F
from hashlib import sha256
from math import lcm
from pathlib import Path

from certify_rotation_inballs import hull
from rational_trig import PI, inv, mul, scale, trig

TAU = F(419, 500)
GRID = 10**8
INNER = F(237, 500)


def body_material(circle, denominator, n):
    vertex_radius = inv(scale(trig(scale(PI, F(1, 2*n)), "cos"), F(2)))
    assert vertex_radius[1] < 1 - INNER
    centers = []
    proxies = []
    for k in range(n):
        angle = scale(PI, F(2*k, n))
        box = tuple(mul(vertex_radius, trig(angle, name)) for name in ("cos", "sin"))
        integer_box = tuple(((lo*GRID).__floor__(), (hi*GRID).__ceil__()) for lo, hi in box)
        centers.append(integer_box)
        proxies.append(tuple((lo + hi)//2 for lo, hi in integer_box))

    def inside(p):
        return all(sum(max(abs(u - lo), abs(u - hi))**2
                       for u, (lo, hi) in zip(p, center)) < GRID**2 for center in centers)

    out = []
    for x, y in circle:
        def point(k):
            return (2*k*x//denominator, 2*k*y//denominator)
        lo, hi = 0, GRID
        while hi-lo > 1:
            mid = (lo+hi)//2
            if inside(point(mid)):
                lo = mid
            else:
                hi = mid
        out.append(point(lo))
    # Arc junctions benefit from adding rational contractions of their proxies.
    out.extend(tuple(u*(GRID-20)//GRID for u in p) for p in proxies)
    assert all(inside(point) for point in out)
    return hull(out), centers


def certify(candidate, sectors, circle, denominator, material, n):
    # Dihedral symmetry reduces directions to [0,pi/n].
    # w in [0,1/3] covers that angular interval.
    angle = scale(PI, F(1, 2*n))
    assert trig(angle, "sin")[1] < trig(angle, "cos")[0]/3
    # The norm of the tan-half-angle derivative is at most 2.
    delta = candidate/(3*sectors)
    contraction = 1-delta/INNER
    worst = None
    grid = contraction.denominator*GRID
    for i in range(sectors):
        w = F(2*i+1, 6*sectors)
        tx = candidate*(1-w*w)/(1+w*w)
        ty = 2*candidate*w/(1+w*w)
        unit = lcm(denominator, grid, tx.denominator, ty.denominator)
        disk_scale = unit//denominator
        body_scale = contraction.numerator*(unit//grid)
        a, b = int(tx*unit), int(ty*unit)
        vertices = hull([(x*disk_scale, y*disk_scale) for x, y in circle]
                        + [(x*body_scale+a, y*body_scale+b) for x, y in material])
        twice = sum(x*v-y*u for (x,y),(u,v) in zip(vertices,vertices[1:]+vertices[:1]))
        area = F(twice, 2*unit**2)
        worst = min(worst, (area, i)) if worst is not None else (area, i)
    return {
        "radius": str(candidate),
        "sectors": sectors,
        "minimum_inner_area": str(worst[0]),
        "minimum_inner_area_display": float(worst[0]),
        "minimum_index": worst[1],
        "verified_exterior_exclusion": worst[0] > TAU,
        "inner_contraction": str(contraction),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, choices=[5,7,9], default=5)
    parser.add_argument("radii", type=F, nargs="+")
    parser.add_argument("--sectors", type=int, default=512)
    args = parser.parse_args()
    raw = Path(__file__).with_name("outer_domain_certificate.json").read_bytes()
    data = json.loads(raw)
    denominator = int(data["disk_polygon_denominator"])
    circle = [tuple(map(int, point)) for point in data["disk_polygon_points"]]
    assert all(4*(x*x+y*y) <= denominator**2 for x,y in circle) and len(circle) == 512
    assert args.sectors>0 and all(0<r<1 for r in args.radii)
    material, centers = body_material(circle, denominator, args.n)
    print(json.dumps({
        "n": args.n,
        "python": platform.python_version(),
        "domain_sha256": sha256(raw).hexdigest(),
        "center_coordinate_enclosures_in_scaled_grid": centers,
        "material_denominator": GRID,
        "material_vertex_count": len(material),
        "tau": str(TAU),
        "entries": [certify(r, args.sectors, circle, denominator, material, args.n) for r in args.radii],
    }, indent=2))


if __name__ == "__main__":
    main()
