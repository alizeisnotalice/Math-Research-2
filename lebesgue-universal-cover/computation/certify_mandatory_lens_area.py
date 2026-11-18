"""Inner triangles and outer tangent triangles for the 13 anchored unit disks."""
from fractions import Fraction
import json

import audit_optimized_eh_chord as g
from certify_eh_general_relocation import rotate_angle
from certify_eh_joint_parameters import specialize


def anchors():
    specialize(Fraction("1.5497"))
    n, m, _, b, _, _, s, t, _ = g.setup()
    r = g.line_pair(m[2], n[1], g.HALF, g.arb(0))
    f = g.unit(5*g.PI/6-g.SIGMA)
    seeds = [s, g.minus(s,n[0]), g.minus(r,f), t,
             (g.arb(".125327173657"), g.arb(".467705236099"))]
    mirrors = [rotate_angle((x,-y),5*g.PI/3) for x,y in seeds]
    return seeds + mirrors + [b, g.minus(b,n[1]), r]


def verify(n):
    disks = [(a, 1-g.dot(a,a)) for a in anchors()]
    assert all(c > 0 for _,c in disks)
    def radial(u):
        values = []
        for center, c in disks:
            projection = g.dot(center,u)
            values.append(projection + (c+projection**2).sqrt())
        return values
    def radius(u):
        vals = radial(u)
        r = vals[0]
        for v in vals[1:]:
            r = r.min(v)
        return r
    step = 2*g.PI/n
    factor = step.sin()/2
    inner, outer = g.arb(0), g.arb(0)
    left = g.unit(g.arb(0))
    first = left_r = radius(left)
    for i in range(n):
        right = g.unit((i+1)*step)
        right_r = first if i+1 == n else radius(right)
        mid = g.unit((2*i+1)*step/2)
        choices = radial(mid)
        j = min(range(len(choices)), key=lambda k: float(choices[k].mid()))
        center = disks[j][0]
        p = tuple(choices[j]*x for x in mid)
        normal = g.minus(p,center)
        height = g.dot(normal,p)
        a, b = g.dot(normal,left), g.dot(normal,right)
        assert a > 0 and b > 0 and height > 0
        inner += factor*left_r*right_r
        outer += factor*height**2/(a*b)
        left, left_r = right, right_r
    assert 0 < inner and inner < outer
    return {"directions":n,"bits":g.ctx.prec,"lower":str(inner.lower()),"upper":str(outer.upper())}


if __name__ == "__main__":
    print(json.dumps(verify(16384), indent=2))
