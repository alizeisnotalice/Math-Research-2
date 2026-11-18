"""Exact perturbation certificate for five curves and an inclusion-minimal eleven-point input."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
from zipfile import ZipFile


def projection(value, axis):
    return sum(a*b for a,b in zip(value, axis))


def distance2(p, q):
    return sum((a-b)**2 for a,b in zip(p,q))


def certify(path):
    base = Path(__file__).resolve().parent
    with ZipFile(base / "global-083732.zip") as source:
        raw = json.loads(source.read("new-points.json"))
        clearance = F(json.loads(source.read("perturbation.json"))["clearance_floor"])
    Q = [[F(x, raw["denominator"]), F(y, raw["denominator"])] for x,y in raw["points"]]
    p, q = Q[4], Q[10]
    axis = tuple(b-a for a,b in zip(p,q))
    normal = (-axis[1], axis[0])
    diameter2 = distance2(p,q)
    assert max(distance2(a,b) for a in Q for b in Q) == diameter2
    assert F(".999999997")**2 < diameter2 < 1
    fixed = []
    for vertex in Q:
        shift = [a-b for a,b in zip(vertex,p)]
        fixed.append((projection(shift,axis)/diameter2, projection(shift,normal)/diameter2))
    assert fixed[4] == (0,0) and fixed[10] == (1,0)
    assert max(distance2(a,b) for a in fixed for b in fixed) == 1

    area = F(".83732")
    target = F(".83730")
    eps = F(".00000298575")
    scale = 1-2*eps
    def area_floor(scale, error):
        return (area*scale**2 - F(22,7)*error**2)/(1+4*error)
    floor = area_floor(scale,eps)
    assert floor > F(".83730000001517") > target
    assert area_floor(1-2*(eps+F(1,10**11)),eps+F(1,10**11)) < target
    # Use an alignment fixing the first diameter endpoint. Every point of Q
    # then moves by at most the expansion of the unit-radius disk.
    normalization = 1/F(".999999997")-1
    assert max(distance2(v,(F(1,2),0)) for v in fixed) <= F(3,4)
    assert normalization + 3*eps < clearance
    assert all(distance2(v,(0,0)) < F(".51")**2 for v in Q)
    assert F("2.02")*eps < clearance
    loose = F(".000005971")
    assert loose < clearance and area_floor(1,loose) > F(".83730000181470")

    centers = [[str(F(1,2)+scale*(x-F(1,2))),str(scale*y)] if i not in (4,10)
               else [str(x),str(y)] for i,(x,y) in enumerate(fixed)]
    ratio_22 = (eps/F(".0000029"))**22
    ratio_18 = (eps/F(".0000027"))**18
    path.write_text(json.dumps({"radius": str(eps), "center_scale": str(scale),
        "fixed_diameter_endpoints": [4,10], "normalized_centers": centers,
        "strict_area_floor": str(floor), "deletion_clearance": str(clearance),
        "maximum_normalized_displacement_upper": str(normalization+3*eps),
        "diameter_squared_before_normalization": str(diameter2),
        "volume_ratios_at_common_target": {"18": str(ratio_18), "22": str(ratio_22)},
        "uncontracted_admissible_radius": str(loose),
        "uncontracted_strict_floor": str(area_floor(1,loose))},indent=2)+"\n")
    print(float(floor), float(ratio_22), float(ratio_18))


if __name__ == "__main__":
    certify(Path(sys.argv[1]))
