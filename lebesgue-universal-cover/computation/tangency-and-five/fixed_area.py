"""Directed area bounds for expanding the four-point cover to contain a unit-diameter regular pentagon."""
import json
from fractions import Fraction as F
from pathlib import Path
import sys
from flint import arb, ctx


def ar(value):
    f=F(value)
    return arb(f.numerator)/f.denominator


def bound(value):
    def rational(side):
        m,e=side.man_exp()
        return str(F(int(m))*2**int(e)) if e>=0 else str(F(int(m),2**(-int(e))))
    return [rational(value.lower()),rational(value.upper())]


def hull(points):
    ordered=sorted(set(points))
    def chain(sequence):
        result=[]
        for x,y in sequence:
            while len(result)>=2:
                a,b=result[-2]; c,d=result[-1]
                if (c-a)*(y-d)>(d-b)*(x-c):break
                result.pop()
            result.append((x,y))
        return result[:-1]
    return chain(ordered)+chain(reversed(ordered))


def certify(path, count=32768):
    ctx.prec=192
    unit=2**56
    pi=arb.pi()
    r=ar("1/2")
    root=arb(3).sqrt()
    def deg(t):return ar(t)*pi/180
    def vec(rad, angle):
        si,co=angle.sin_cos()
        return [rad*co,rad*si]
    radius=1/(2*deg(18).cos())
    b=radius*(root*deg(6).cos()/2+deg(72).cos())
    d=(radius*deg(6).sin())**2
    t=d/(b+(b*b-d/4).sqrt())
    delta=radius*deg(6).cos()+root*t/2-r
    area0=pi/8+root/4
    cap=r*delta**2/((r+delta)**2-r**2).sqrt()
    assert cap>0
    p=vec(t,deg(330))
    def integer(v):
        x,y=v
        assert x*x+y*y<ar(".60")**2
        def lattice(a):
            k=(a*unit+ar(".5")).floor().unique_fmpz()
            assert k is not None and abs(a*unit-int(k))<ar(".5")
            return int(k)
        return (lattice(x),lattice(y))
    circle=[integer(vec(r,2*pi*j/count)) for j in range(count)]
    triangle=[integer(vec(1/root,deg(30+120*j))) for j in range(3)]
    pentagon=[integer([x+y for x,y in zip(p,vec(radius,deg(6+72*j)))]) for j in range(5)]
    # A second feasible placement uses the pentagon locked inside the Gibbs two-cap host.
    lock_shift=2*(radius*deg(6).cos()-r)/root
    locked=[integer([v[0],v[1]+lock_shift]) for j in range(5)
            for v in [vec(radius,deg(18+72*j))]]
    def poly_area(vertices):
        return ar(F(sum(a*d-b*c for (a,b),(c,d) in zip(vertices,vertices[1:]+vertices[:1])),2*unit**2))
    polygon=hull(circle+triangle+pentagon)
    base=poly_area(polygon)
    lock_area=poly_area(hull(circle+triangle+locked))
    # Both true and rounded convex hulls are in radius .61 and have perimeter <4.
    error=(1-(pi/count).cos())/2+ar(1)/unit
    assert 2*pi*ar(".61")<4
    shift=4*error+pi*error**2
    fixed=(base-shift).union(base+shift)
    lock_interval=(lock_area-shift).union(lock_area+shift)
    assert area0+cap < lock_interval < fixed
    out=dict(precision_bits=ctx.prec, circle_sides=count,lattice_denominator=unit,
             fixed_hull_vertices=len(polygon),base_area=bound(area0),disk_tangent_cap=bound(cap),
             all_poses_lower=bound(area0+cap),distance_minimizing_pose_area=bound(fixed), locked_pentagon_pose_area=bound(lock_interval),
             min_area_gain_lower=bound(cap),min_area_gain_upper=bound(lock_interval-area0))
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(float(area0),float(cap),float((area0+cap).lower()),float(fixed.lower()),float(fixed.upper()),float(lock_interval.lower()),float(lock_interval.upper()))


if __name__=="__main__":
    certify(Path(sys.argv[1]))
