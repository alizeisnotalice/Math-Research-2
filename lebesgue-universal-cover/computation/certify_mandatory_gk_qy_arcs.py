"""Certify two symmetric nine-anchor lenses, two full Sprague arcs and the AH circular arc."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import audit_optimized_eh_chord as g
from audit_eh_continuous_placement import project
from certify_eh_general_relocation import exact_fraction,rotate_angle
from certify_eh_joint_parameters import specialize


def verify(name):
    specialize(Fraction('1.5497'))
    n,m,*_=g.setup()
    if name=='GK':
        s=1
        c=g.line_pair(n[3],m[2])
        l=g.circle_line(c,n[5],1)
        r=project(c,m[5])
        pts={'C':c,'E':g.plus(c,n[0]),'G':l,'K':r,'A':g.scale(g.arb('.5634'),g.unit(g.PI/6))}
        reflected=('C','E','G','K')
        unit=[('C','E'),('C','G'),('C','K'),('G','JG'),('JC','JE'),('JC','JG'),('JC','JK')]
        start=0
        three=('E','JG','JC')
        strict=[(0,'A'),(1,'JG')]
        equality=(2,'C')
    else:
        s=3
        c=g.line_pair(n[4],m[3])
        l=project(c,n[1])
        r=(g.arb(0),c[1]+(1-c[0]**2).sqrt())
        pts={'C':c,'Q':l,'L':g.plus(c,m[0]),'Y':r,'A':(g.HALF,g.arb('.038'))}
        reflected=('C','Q','L','A')
        unit=[('C','Q'),('C','L'),('C','Y'),('Y','JC'),('A','JA'),('JC','JQ'),('JC','JL')]
        start=1
        three=('JA','JC','Q')
        strict=[(1,'Y'),(2,'JL')]
        equality=(0,'L')
    mirror=lambda z:rotate_angle((z[0],-z[1]),s*g.PI/3)
    pts.update({'J'+k:mirror(pts[k]) for k in reflected})
    old={frozenset(v) for v in unit}
    gaps=[(1-g.dot(g.minus(p,q),g.minus(p,q))).lower() for (i,p),(j,q) in combinations(pts.items(),2)
          if frozenset((i,j)) not in old]
    assert min(gaps)>0
    lam=g.arb('.0000001')
    def angles(t):
        return [start*g.PI+t+i*2*g.PI/3 for i in range(3)]
    rad=lambda d:exact_fraction(d)*g.PI/180
    # Contacts in 'three' have projection exactly 1/2 at the zero angle.
    b=sum(g.dot(g.rot(g.unit(t)),pts[k]) for t,k in zip(angles(g.arb(0)),three,strict=True))
    angle=rad(Fraction(1))
    deriv=(1-lam)*(b*angle.cos()-3*g.HALF*angle.sin())-9*lam/100
    assert deriv>0
    cells=[]
    def scan(a,b):
        window=rad(a).union(rad(b))
        excess=sum(max(g.dot(g.unit(t),p).lower() for p in pts.values()) for t in angles(window))-3*g.HALF
        margin=(1-lam)*excess-3*lam/100
        if margin>0:
            cells.append((str(a),str(b),str(margin.lower())))
            return
        assert b-a>Fraction(1,4096)
        mid=(a+b)/2
        scan(a,mid)
        scan(mid,b)
    scan(Fraction(1),Fraction(30))
    assert min((1-lam)*(g.dot(m[i],pts[k])-g.HALF)-lam/100 for i,k in strict)>0
    edge_coefficient=(3*g.SIGMA).cos()/100
    assert edge_coefficient>0 and equality[0]%2==0
    positive={0,1,2}
    prime={(s-k-1)%6 for k in positive}
    poses=[]
    for j in range(6):
        if not ({(3-j)%6,(5-j)%6}&positive):poses.append([j,False])
        if not ({(j-4)%6,(j-6)%6}&prime):poses.append([j,True])
    assert poses==[[0,False],[s,True]]
    left,right=g.minus(l,c),g.minus(r,c)
    assert g.cross(left,right)>0 and g.dot(left,right)>0
    span=g.arb.atan2(g.cross(left,right),g.dot(left,right))
    out={'anchor_coordinates':{k:list(map(str,z)) for k,z in pts.items()},'exact_unit_pairs':unit,
         'other_pair_count':len(gaps),'other_squared_diameter_slack_lower':str(min(gaps)),
         'first_degree_derivative_lower':str(deriv.lower()),'angle_cells':cells,
         'interior_slack_lower':str(min(g.arb(e[2]).lower() for e in cells)),'allowed_poses':poses,
         'strict_cap_margins':{str(i):str(((1-lam)*(g.dot(m[i],pts[k])-g.HALF)-lam/100).lower()) for i,k in strict},
         'equality_cap':equality,'equality_cap_excess_per_mixing_weight':str(edge_coefficient.lower()),
         'whole_arc_angle_radians':str(span)}
    if name=='QY':
        center=mirror(c)
        midpoint=(-g.HALF/2,-g.arb(3).sqrt()/4)
        rr=g.arb(3).sqrt()
        l0=(rr/2-g.HALF,(rr/2-g.HALF-1)/rr)
        z=g.circle_pair(l0,midpoint,-1)
        v=g.minus(z,center)
        vv=g.dot(v,v)
        ta=g.plus(center,g.scale(1/vv,g.minus(v,g.scale((vv-1).sqrt(),g.rot(v)))))
        u0=g.minus(r,center)
        u1=g.minus(ta,center)
        u2=g.minus(mirror(l),center)
        turns=[g.cross(u0,u1),g.cross(u1,u2),g.dot(u0,u2)]
        assert min(turns)>0
        out['AH_Y_Ta_sweep_radians']=str(g.arb.atan2(turns[0],g.dot(u0,u1)))
        out['AH_radial_order_slacks']=list(map(str,turns))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    out={'sigma_degrees':'1.5497','precision_bits':g.ctx.prec,'cases':{tag:verify(tag) for tag in ('GK','QY')}}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
