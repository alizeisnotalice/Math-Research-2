"""Certify positive-weight forcing and a protrusion beyond the old disk/lens hull."""
import json
from fractions import Fraction as F

from certify_mandatory_lens_area import anchors, g

WEIGHT = F(147, 2500)


def verify():
    points = anchors()
    lam = g.exact(WEIGHT)
    def angle(d):
        return g.exact(d)*g.PI/180

    def normals(t):
        return [g.unit(t+g.PI+2*i*g.PI/3) for i in range(3)]

    # E, QS, B achieve the three half-widths at zero.
    velocity = sum(g.dot(g.rot(u),points[k]) for u,k in zip(normals(g.arb(0)),(1,5,10)))
    derivative = (1-lam)*(velocity*angle(F(1)).cos()-3*g.HALF*angle(F(1)).sin())-9*lam/100
    assert derivative > 0
    bounds = []
    def cover(a,b):
        t = angle(a).union(angle(b))
        lower = sum(max(g.dot(u,z).lower() for z in points) for u in normals(t))-3*g.HALF
        mixed = (1-lam)*lower-3*lam*(3*t).sin()/100
        if mixed > 0:
            bounds.append(mixed.lower())
            return
        assert b-a > F(1,2**23)
        m = (a+b)/2
        cover(a,m)
        cover(m,b)
    cover(F(1),F(30))
    _, m, _, _, _, _, _, _, _ = g.setup()
    cap = (1-lam)*(g.dot(m[1],points[8])-g.HALF)-lam/100
    assert cap > 0 and (3*g.SIGMA).cos() > 0

    # For L in B(U,1), n=R-U is the unit outward normal at R.
    r, u = points[12], points[2]
    n = g.minus(r,u)
    hx = g.dot(n,r)
    third_sine = 3*n[1]-4*n[1]**3
    hh = g.HALF+third_sine/100
    old = hx.max(g.HALF)
    excess = (1-lam)*hx+lam*hh-old
    assert excess > g.arb(".00058606535")
    # The old hull is in x <= old. A base chord of the required disk at
    # x=.4995 and the new point make a triangle with a certified outer tip.
    x0 = g.arb(".4995")
    base_half = (g.HALF**2-x0**2).sqrt()
    area = base_half * excess**2/(old+excess-x0)
    assert area > g.arb(".00000706988")
    return {"weight":str(WEIGHT),"bits":g.ctx.prec,"first_degree_derivative":str(derivative),
            "root_cells":len(bounds),"minimum_root_lower":str(min(bounds)),
            "middle_cap_margin":str(cap),"excess_over_disk_and_lens":str(excess),
            "additional_area_lower":str(area.lower())}


if __name__ == "__main__":
    print(json.dumps(verify(),indent=2))
