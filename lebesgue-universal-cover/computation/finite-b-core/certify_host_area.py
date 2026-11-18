"""Enclose the convex hull of K0 and 150 specifically placed whole width-one bodies."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

import flint
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from certify_nonmonotone_ah_core import double_area, hull

TARGET = Fraction(".844082046521218581802948589770767704354")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def endpoints(x):
    def dyadic(y):
        m, e = map(int, y.man_exp())
        return Fraction(m) * Fraction(2) ** e
    return dyadic(x.lower()), dyadic(x.upper())


def prove(poses_path, p_path, c_path, core_path):
    ctx.prec = 192
    paths = (poses_path, p_path, c_path, core_path)
    poses, pentagons, curves, core = [json.loads(p.read_bytes()) for p in paths]
    assert pentagons["poses_sha256"] == curves["poses_sha256"] == digest(poses_path)
    assert pentagons["script_sha256"] == digest(HERE / "sample_pentagons.py")
    assert curves["source_sha256"] == digest(HERE / "sample_profiles.py")
    for deps in (pentagons["source_sha256"], curves["input_sha256"], curves["dependency_sha256"]):
        for name, checksum in deps.items():
            assert digest(HERE.parent / name) == checksum
    records = pentagons["per_body"] + curves["bodies"]
    assert len(records) == 150 and len(pentagons["per_body"]) == 117
    for actual, expected in zip(sorted(records, key=lambda r: r["index"]), poses["poses"], strict=True):
        assert all(actual[k] == expected[k] for k in ("index","source","body_index","turn","reflection"))
    assert sorted(r["index"] for r in records) == list(range(150))
    n = pentagons["directions"]
    assert n == curves["N"] and n >= 4 and curves["complete"]
    assert all(r.get("samples", r.get("count")) == n for r in records)
    assert pentagons["sample_count"] + curves["sample_count"] == 150 * n
    grid = core["vertex_grid"]
    scale, rem = divmod(grid, curves["scale"])
    assert scale > 0 and rem == 0 and curves["scale"] == pentagons["grid"]
    coordinate_errors = (Fraction(pentagons["coordinate_rounding_error_strict"]),
                         Fraction(curves["coordinate_error_asserted"]))
    grid_error = Fraction(3, 10**11)
    assert 2 * max(coordinate_errors)**2 < grid_error**2

    core_vertices = list(map(tuple, core["outer_integer_vertices"]))
    assert double_area(core_vertices) == int(core["double_area_integer"])
    points = core_vertices.copy()
    for x, y in pentagons["vertices_ccw_integer"] + curves["combined_hull"]:
        points.append((scale*x, scale*y))
    polygon = hull(points)
    area = Fraction(double_area(polygon), 2*grid**2)
    perimeter = sum((arb((x-u)**2 + (y-v)**2).sqrt()
                     for (x,y),(u,v) in zip(polygon, polygon[1:]+polygon[:1])), arb(0))/grid
    r = 2*(arb.pi()/(2*n)).sin()**2 + arb(grid_error.numerator)/grid_error.denominator
    added = r*perimeter + arb.pi()*r*r
    upper = area + endpoints(added)[1]
    assert upper < TARGET
    # Rational endpoints readable directly, retaining the direction of rounding.
    decimal_scale = 10**35
    ceiling = (-(-upper.numerator*decimal_scale//upper.denominator))
    gap = TARGET - upper
    return {
        "theorem": "area of the convex hull of K0 and all 150 chosen copies",
        "input_sha256": {p.name: digest(p) for p in paths},
        "python": platform.python_version(), "flint": flint.__version__, "bits": ctx.prec,
        "body_count": len(records), "N": n, "sample_count": 150*n,
        "outer_core_vertices": len(core_vertices), "integer_vertices": len(polygon),
        "integer_vertex_sha256": sha256("".join(f"{x},{y}\n" for x,y in polygon).encode()).hexdigest(),
        "integer_grid": grid, "polygon_area": str(area),
        "polygon_perimeter": list(map(str, endpoints(perimeter))),
        "point_error_upper": str(grid_error), "enlargement_radius": list(map(str, endpoints(r))),
        "steiner_correction": list(map(str, endpoints(added))),
        "host_area_upper": str(upper), "host_area_decimal_upper": f"{ceiling//decimal_scale}.{ceiling%decimal_scale:035d}",
        "B": str(TARGET), "slack_lower": str(gap),
    }


if __name__ == "__main__":
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("poses", "pentagons", "curves", "core", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    result = prove(args.poses, args.pentagons, args.curves, args.core)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(result["host_area_decimal_upper"])
