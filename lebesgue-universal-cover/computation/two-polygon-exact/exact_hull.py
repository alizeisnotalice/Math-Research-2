"""Check the quadratic area from the exact convex hull of all 41 placed points."""
from fractions import Fraction as F
from functools import cmp_to_key
import json
from math import isqrt
from pathlib import Path
import sys


def verify(work):
    data = json.loads((work / "micro-scan.json").read_bytes())
    value = json.loads((work / "algebraic.json").read_bytes())
    a, b, c = map(F, value["u_quadratic_coefficients"])
    disc = b*b - 4*a*c
    assert disc > 0 and F(value["b"]) != 0
    assert isqrt(disc.numerator)**2 != disc.numerator or isqrt(disc.denominator)**2 != disc.denominator
    lo, hi = F("0.10859334269301195844"), F("0.10859334269301195845")
    poly = lambda x: (a*x+b)*x+c
    assert poly(lo)*poly(hi) < 0 and (2*a*lo+b)*(2*a*hi+b) > 0
    u2, u1 = -c/a, -b/a

    def add(x, y):
        return x[0]+y[0], x[1]+y[1]

    def sub(x, y):
        return x[0]-y[0], x[1]-y[1]

    def mul(x, y):
        return x[0]*y[0]+x[1]*y[1]*u2, x[0]*y[1]+x[1]*y[0]+x[1]*y[1]*u1

    def divide(x, y):
        conjugate = (y[0]+y[1]*u1, -y[1])
        norm, zero = mul(y, conjugate)
        assert zero == 0 and norm
        n = mul(x, conjugate)
        return n[0]/norm, n[1]/norm

    def sign(x):
        if x == (0, 0):
            return 0
        low, high = lo, hi
        while True:
            left, right = x[0]+x[1]*low, x[0]+x[1]*high
            if left > 0 and right > 0:
                return 1
            if left < 0 and right < 0:
                return -1
            mid = (low+high)/2
            if poly(mid)*poly(low) < 0:
                high = mid
            else:
                low = mid

    def cross(x, y):
        return sub(mul(x[0], y[1]), mul(x[1], y[0]))

    def vector(a, b):
        return sub(a[0], b[0]), sub(a[1], b[1])

    def dot(x, y):
        return add(mul(x[0], y[0]), mul(x[1], y[1]))

    unit, u = (F(1), F(0)), (F(0), F(1))
    sq = mul(u, u)
    co = divide(sub(unit, sq), add(unit, sq))
    si = divide(add(u, u), add(unit, sq))
    p = [((F(x), 0), (F(y), 0)) for x, y in data["hullP"]]
    q = [(sub(mul((F(x), 0), co), mul((F(y), 0), si)),
          add(mul((F(x), 0), si), mul((F(y), 0), co))) for x, y in data["hullQ"]]
    normals, offsets = [], []
    for edge, vertex in ((0, 5), (7, 8), (26, 3)):
        dx, dy = vector(p[(edge+1) % len(p)], p[edge])
        n = dy, (-dx[0], -dx[1])
        normals.append(n)
        offsets.append(dot(n, vector(p[edge], q[vertex])))
    n, m, k = normals
    det = cross(n, m)
    tx = divide(sub(mul(offsets[0], m[1]), mul(n[1], offsets[1])), det)
    ty = divide(sub(mul(n[0], offsets[1]), mul(offsets[0], m[0])), det)
    assert dot(k, (tx, ty)) == offsets[2]
    points = p + [(add(x, tx), add(y, ty)) for x, y in q]
    ordered = sorted(points, key=cmp_to_key(lambda x,y: sign(sub(x[0],y[0])) or sign(sub(x[1],y[1]))))

    def chain(seq):
        path = []
        for z in seq:
            while len(path) >= 2 and sign(cross(vector(path[-1], path[-2]), vector(z, path[-1]))) <= 0:
                path.pop()
            path.append(z)
        return path[:-1]

    boundary = chain(ordered) + chain(reversed(ordered))
    area = (F(0), F(0))
    for x, y in zip(boundary, boundary[1:]+boundary[:1]):
        area = add(area, cross(x, y))
    assert (area[0]/2, area[1]/2) == (F(value["a"]), F(value["b"]))
    return {"input_points": len(points), "hull_vertices": len(boundary),
            "quadratic_degree": 2, "exact_shoelace_agrees": True}


if __name__ == "__main__":
    print(json.dumps(verify(Path(sys.argv[1])), indent=2))
