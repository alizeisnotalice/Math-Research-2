#!/usr/bin/env python3
"""Certify a complete obstruction to fitting a width-one body in a fixed EH cover."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_all_completions import Lens
from certify_eh_general_relocation import exact_fraction,rotate_angle


def verify(case):
    g.SIGMA=exact_fraction(Fraction(case['sigma']))*g.PI/180
    g.trial(case['delta'],tight=True)
    delta=g.arb(case['delta'])
    phi=exact_fraction(Fraction(case['support_angle_degrees']))*g.PI/180
    n,m,c3,b3,_,u,_,t,j=g.setup()
    axis=n[5]

    def reflect(p):
        return g.minus(g.scale(2*g.dot(axis,p),axis),p)

    def sq(p):
        return g.dot(p,p)

    raw=list(map(g.arb,case['coordinates']))
    z,zp=tuple(raw[:2]),tuple(raw[2:4])
    b,bp=(g.plus(c3,(g.arb(0),height)) for height in raw[4:6])
    w,wp=tuple(raw[6:8]),tuple(raw[8:10])
    q=g.plus(g.scale(g.HALF,axis),g.scale(raw[10],g.rot(axis)))
    points={'z':z,'b':b,'b_op':g.minus(b,n[0]),'w':w,'w_op':g.plus(w,m[5]),
            'q':q,'q_op':g.minus(q,axis),'z_Q':reflect(zp),'b_Q':reflect(bp),
            'b_Q_op':reflect(g.minus(bp,n[0])),'w_Q':reflect(wp),'w_Q_op':reflect(g.plus(wp,m[5]))}
    unit_pairs=[('b','b_op'),('w','w_op'),('q','q_op'),('b_Q','b_Q_op'),('w_Q','w_Q_op')]
    for j0,i in enumerate((11,13,15)):
        contact=tuple(raw[i:i+2]); direction=g.scale((-1)**j0,g.unit(phi+j0*g.PI/3))
        name=f'contact_{j0}'
        points[name]=contact
        points[name+'_op']=g.minus(contact,direction)
        unit_pairs.append((name,name+'_op'))

    exact={frozenset(x) for x in unit_pairs}
    gaps=[1-sq(g.minus(p,v)) for (a,p),(b,v) in combinations(points.items(),2) if frozenset((a,b)) not in exact]
    assert all(slack>0 for slack in gaps)
    anchors={'m5_opposite':g.dot(g.neg(m[5]),w)-g.HALF,
             'Qm5_opposite':g.dot(g.neg(m[5]),wp)-g.HALF,
             'm3_opposite':g.dot(m[0],points['z_Q'])-g.HALF,
             'Qm3_opposite':g.dot(m[0],reflect(z))-g.HALF}
    assert all(v>0 for v in anchors.values())

    co,si=g.unit(delta)
    radial=g.minus(u,c3)
    u1=g.plus(c3,g.minus(g.scale(co,radial),g.scale(si,g.rot(radial))))
    chord=g.minus(t,u1)
    assert sq(g.minus(u1,b3))<1 and sq(g.minus(t,c3))<1
    assert g.cross(chord,g.minus(j,u1))<0
    assert g.cross(chord,g.minus(g.circle_pair(b3,c3,1),u1))>0
    assert g.dot(chord,g.minus(u1,c3))<0 and g.dot(chord,g.minus(t,b3))>0
    ear=[]
    for z0 in (z,zp):
        ear.extend((1-sq(g.minus(z0,c3)),1-sq(g.minus(z0,b3)),-g.cross(chord,g.minus(z0,u1))))
    assert all(v>0 for v in ear)

    # Unit contacts give F'(0) exactly for every completion. Its derivative
    # varies by at most 7.8 times the angle variation in radians.
    derivative=(g.dot(g.rot(n[0]),b)-g.dot(g.rot(n[1]),points['b_Q_op'])
                +g.dot(g.rot(n[2]),points['q_op']))
    epsilon=Fraction(1,4)
    lower_derivative=derivative-g.arb('7.8')*exact_fraction(epsilon)*g.PI/180
    assert lower_derivative>0

    lens=Lens(points)
    cache={}

    def negative_support_sum(a):
        if a not in cache:
            angle=exact_fraction(a)*g.PI/180
            normals=[g.scale((-1)**(k+1),g.unit(angle+k*g.PI/3)) for k in range(3)]
            cache[a]=sum(lens.upper_support(v) for v in normals).upper()
        return cache[a]

    cells=[]

    def partition(a,b,depth):
        h=exact_fraction(b-a)*g.PI/360
        # Endpoint normals use the same nonnegative angular interpolation.
        # The endpoint upper bounds are positive, so use the secant factor.
        upper=max(negative_support_sum(a),negative_support_sum(b))
        assert upper>0 and h.cos()>0
        lower=g.arb('1.5')-upper/h.cos()
        if lower>0:
            cells.append((str(a),str(b),lower.lower()))
            return
        assert depth<14,(str(a),str(b),str(lower))
        middle=(a+b)/2
        partition(a,middle,depth+1)
        partition(middle,b,depth+1)

    partition(epsilon,60-epsilon,0)
    excluded={}
    for shift in range(6):
        for reflection in (False,True):
            if (shift,reflection) in ((0,False),(4,True)):
                continue
            claims=[]
            for name,(x,y) in points.items():
                image=rotate_angle((x,-y if reflection else y),shift*g.PI/3)
                for cap in (3,5):
                    claims.append((g.dot(m[cap],image)-g.HALF,name,cap))
            margin,name,cap=max(claims,key=lambda item:item[0].lower())
            assert margin>0,(shift,reflection,str(margin))
            excluded[f'{shift}_{int(reflection)}']=[name,cap,str(margin.lower())]
    least,angle_margin=min(gaps,key=lambda z:z.lower()),min(r[2] for r in cells)
    return {'sigma':case['sigma'],'delta':case['delta'],'num_witness_points':len(points),
            'exact_unit_pairs':unit_pairs,'squared_slack_other_pairs':str(least.lower()),
            'strict_ear_slack':str(min(ear,key=lambda z:z.lower()).lower()),
            'uniform_F_prime_lower_first_quarter_degree':str(lower_derivative.lower()),
            'uniform_F_lower_interior':str(angle_margin),'interior_angle_cell_count':len(cells),
            'F_positive_sector_degrees':['0','60'],'only_d6_poses_in_P':['identity','Q2'],
            'excluded_d6_poses':excluded,'angle_cells':[[a,b,str(v)] for a,b,v in cells]}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    cases=json.loads(args.input.read_bytes())
    reports=[verify(case) for case in cases]
    args.output.write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in row.items() if k not in ('angle_cells','excluded_d6_poses','exact_unit_pairs')} for row in reports],indent=2))
