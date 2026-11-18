"""Exact polygon area bounds for the convex hull of the three forced lenses and disk."""
import json
from fractions import Fraction as F

from certify_mandatory_lens_area import anchors, g
from certify_mandatory_gk_qy_arcs import verify as arc_verify
from certify_nine import S, append_point, hull


def twice_area(points):
    poly = hull(points)
    return sum(x*v-y*u for (x,y),(u,v) in zip(poly,poly[1:]+poly[:1]))


def verify(n):
    groups = [anchors()]
    for name in ("GK","QY"):
        data = arc_verify(name)
        groups.append([tuple(map(g.arb,p)) for p in data["anchor_coordinates"].values()])
    disks = [[(p,1-g.dot(p,p)) for p in row] for row in groups]
    assert all(c > 0 for row in disks for _,c in row)
    inner, outer = [], [(x*S//4,y*S//4) for x in (-1,1) for y in (-1,1)]
    step = 2*g.PI/n
    for i in range(n):
        left,right,mid = [g.unit(step*k) for k in (i,i+1,g.arb(i)+g.HALF)]
        def add(points,u,r):
            append_point(points,u[0]*r,u[1]*r)
        add(inner,left,g.HALF)
        add(outer,left,g.HALF/g.dot(mid,left))
        add(outer,right,g.HALF/g.dot(mid,right))
        for row in disks:
            def radial(u):
                proj = [g.dot(p,u) for p,_ in row]
                return [d+(d*d+c).sqrt() for d,(_,c) in zip(proj,row)]
            vals = radial(left)
            r = vals[0]
            for v in vals[1:]:
                r = r.min(v)
            add(inner,left,r)
            rm = radial(mid)
            j = min(range(len(rm)),key=lambda k:float(rm[k].mid()))
            contact = g.scale(rm[j],mid)
            normal = g.minus(contact,row[j][0])
            h = g.dot(normal,contact)
            for u in (left,right):
                denom = g.dot(normal,u)
                assert denom > 0 and h > 0
                add(outer,u,h/denom)
    low = F(twice_area(inner),2*(S+4)**2)
    high = F(twice_area(outer)*(S+8)**2,2*S**4)
    assert 0 < low < high < 1
    return {"sectors":n,"lower":str(low),"upper":str(high),"display":[float(low),float(high)]}


if __name__ == "__main__":
    print(json.dumps(verify(4096), indent=2))
