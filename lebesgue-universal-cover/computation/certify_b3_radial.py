"""Exact exterior exclusion for the entire B3 translation disk via parallel chords."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
import platform

from certify_rotation_inballs import HERE, TAU, hull

GRID = 10**8
RADIUS = F(7,100)
SECTORS = 1024
INNER_RADIUS = F(21,50)


def material(den,circle):
    # Each true center coordinate lies within one unit of the integer proxy /GRID.
    # floor(GRID/sqrt(3)) and floor(GRID/(2sqrt(3))) via integer squares.
    c = isqrt(GRID*GRID//3)
    assert 3*c*c <= GRID**2 < 3*(c+1)**2
    centers=[(c,0),(-c//2,GRID//2),(-c//2,-GRID//2)]
    def belongs(p):
        return all((abs(p[0]-x)+1)**2+(abs(p[1]-y)+1)**2<=GRID**2 for x,y in centers)
    out=[]
    for x,y in circle:
        def point(k):
            return (2*k*x//den,2*k*y//den)
        # Bisection only chooses material; the final test verifies membership directly.
        lo,hi=0,GRID
        while hi-lo>1:
            mid=(lo+hi)//2
            if belongs(point(mid)):
                lo=mid
            else:
                hi=mid
        z=point(lo)
        assert belongs(z)
        out.append(z)
    # Three slightly contracted corners repair irrational vertex-direction sampling.
    for x,y in centers:
        z=(x*(GRID-20)//GRID,y*(GRID-20)//GRID)
        assert belongs(z)
        out.append(z)
    return hull(out)


def run():
    domain = HERE/'outer_domain_certificate.json'
    raw=domain.read_bytes(); source=json.loads(raw)
    den=int(source['disk_polygon_denominator']); circle=[tuple(map(int,p)) for p in source['disk_polygon_points']]
    assert len(circle)==512 and all(4*(x*x+y*y)<=den**2 for x,y in circle)
    body=material(den,circle)
    assert INNER_RADIUS < 1-F(57735028,10**8)
    assert 3*57735028**2>10**16
    # tan-half parametrization over 0..1 sweeps the first quadrant. Derivative norm <= 2.
    # Every angle (modulo dihedral B3 symmetry) has a representative in that quadrant.
    delta=RADIUS/SECTORS
    contraction=1-delta/INNER_RADIUS
    assert 0<contraction<1
    worst=None
    for i in range(SECTORS):
        w=F(2*i+1,2*SECTORS)
        tx=RADIUS*(1-w*w)/(1+w*w); ty=RADIUS*2*w/(1+w*w)
        # All coordinates on one exact integer grid per sector.
        from math import lcm
        unit=lcm(den,contraction.denominator*GRID,tx.denominator,ty.denominator)
        disk=[(x*(unit//den),y*(unit//den)) for x,y in circle]
        k=contraction.numerator*(unit//(contraction.denominator*GRID))
        a=tx.numerator*(unit//tx.denominator);b=ty.numerator*(unit//ty.denominator)
        p=hull(disk+[(k*x+a,k*y+b) for x,y in body]); twice=sum(x*v-y*u for (x,y),(u,v) in zip(p,p[1:]+p[:1]))
        area=F(twice,2*unit*unit)
        assert area>TAU
        if worst is None or area<worst[0]:
            worst=(area,i)
    return {'python':platform.python_version(),'domain_sha256':sha256(raw).hexdigest(),'sectors':SECTORS,
            'material_vertices':len(body),'radius':str(RADIUS),'w_range':[0,1],
            'delta':str(delta),'body_contraction':str(contraction),'threshold':str(TAU),
            'minimum_inner_area':str(worst[0]),'minimum_display':float(worst[0]),'minimizing_sector_index':worst[1],
            'scope':'All normalized F8 configurations of area <=419/500 must have |t_B3|<7/100.'}

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
