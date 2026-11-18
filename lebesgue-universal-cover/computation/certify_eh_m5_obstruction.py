#!/usr/bin/env python3
"""Certify finite witnesses forcing failure of the R90 right-root pose on m5."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import audit_optimized_eh_chord as geo
from certify_eh_all_completions import Lens
from certify_eh_general_relocation import exact_fraction,rotate_angle


def candidate(case):
    sigma=Fraction(case['sigma'])
    geo.SIGMA=exact_fraction(sigma)*geo.PI/180
    delta,phi=map(geo.arb,(case['delta'],case['phi']))
    geo.trial(case['delta'],tight=True)
    n,m,c3,b3,_,u,_,t,j=geo.setup()
    arb=geo.arb
    dot,plus,minus,neg,scale,rot=(geo.dot,geo.plus,geo.minus,geo.neg,geo.scale,geo.rot)
    axis=n[5]

    def reflect(y):
        return minus(scale(2*dot(axis,y),axis),y)

    theta=(arb(240)-phi)*geo.PI/180
    v=[rotate_angle(w,-theta) for w in n[:3]]
    raw=[arb(value) for value in case['coordinates']]
    z,zp=tuple(raw[:2]),tuple(raw[2:4])
    b,bp=plus(c3,(arb(0),raw[4])),plus(c3,(arb(0),raw[5]))
    w,wp=tuple(raw[6:8]),tuple(raw[8:10])
    q=plus(scale(geo.HALF,axis),scale(raw[10],rot(axis)))
    p0,p1,p2,x=(tuple(raw[i:i+2]) for i in (11,13,15,17))
    # Eliminate the compatibility equation exactly by changing only
    # the v2 component of p2. The three vj are exact unit vectors.
    correction=geo.HALF-dot(v[0],p0)+dot(v[1],p1)-dot(v[2],p2)
    p2=plus(p2,scale(correction,v[2]))
    P={'z':z,'b':b,'b_opposite':minus(b,n[0]),'w':w,'w_opposite':plus(w,m[5]),
       'q':q,'q_opposite':minus(q,axis),'z_Q':reflect(zp),'b_Q':reflect(bp),
       'b_Q_opposite':reflect(minus(bp,n[0])),'w_Q':reflect(wp),
       'w_Q_opposite':reflect(plus(wp,m[5])),'p0':p0,'p1':p1,'p2':p2,'x':x,
       'p0_opposite':minus(p0,v[0]),'p1_opposite':minus(p1,v[1]),'p2_opposite':minus(p2,v[2])}
    equal=[('b','b_opposite'),('w','w_opposite'),('q','q_opposite'),('b_Q','b_Q_opposite'),
           ('w_Q','w_Q_opposite'),('p0','p0_opposite'),('p1','p1_opposite'),('p2','p2_opposite')]

    def squared(y):
        return dot(y,y)

    unit_pairs={frozenset(pair) for pair in equal}
    slacks=[1-squared(minus(a,b)) for (i,a),(j,b) in combinations(P.items(),2)
            if frozenset((i,j)) not in unit_pairs]
    assert all(value>0 for value in slacks)

    # These opposite-strip anchors keep every constant-width completion in P and QP.
    anchors={
        'm5_opposite':dot(neg(m[5]),w)-geo.HALF,
        'Qm5_opposite':dot(neg(m[5]),wp)-geo.HALF,
        'm3_opposite':dot(m[0],P['z_Q'])-geo.HALF,
        'Qm3_opposite':dot(m[0],reflect(z))-geo.HALF,
    }
    assert all(value>0 for value in anchors.values())

    co,si=geo.unit(delta)
    radial=minus(u,c3)
    u1=plus(c3,minus(scale(co,radial),scale(si,rot(radial))))
    chord=minus(t,u1)
    assert squared(minus(u1,b3))<1 and squared(minus(t,c3))<1
    assert geo.cross(chord,minus(j,u1))<0
    assert geo.cross(chord,minus(geo.circle_pair(b3,c3,1),u1))>0
    assert dot(chord,minus(u1,c3))<0 and dot(chord,minus(t,b3))>0
    ear_slacks=[]
    for hit in (z,zp):
        ear_slacks.extend((1-squared(minus(hit,c3)),1-squared(minus(hit,b3)),
                           -geo.cross(chord,minus(hit,u1))))
    assert all(value>0 for value in ear_slacks)

    a0,a1=geo.HALF-dot(v[0],p0),geo.HALF-dot(v[1],p1)
    translation=(a0,(2*a1-a0)/arb(3).sqrt())
    violation=dot(m[5],plus(rotate_angle(x,theta),translation))-geo.HALF
    assert violation>0

    # At the root F has period 120 degrees. Estimate its derivative near
    # this particular exact triple of support contacts using their unit mates.
    derivative=sum(sign*dot(rot(u),p) for sign,u,p in zip((1,-1,1),v,(p0,p1,p2)))
    neighborhood=arb('0.3')*geo.PI/180
    derivative_upper=derivative+arb('7.8')*neighborhood
    assert derivative_upper<0

    lens=Lens(P)
    left=arb(30)
    end=phi-arb('0.3')
    width=(end-left)/128
    assert width>0
    signs=[]
    for i in range(128):
        angle=(left+i*width).union(left+(i+1)*width)*geo.PI/180
        vec=[geo.unit(angle+k*geo.PI/3) for k in range(3)]
        lower=arb('1.5')-sum(lens.upper_support(scale(sign,d)) for sign,d in zip((-1,1,-1),vec))
        signs.append(lower)
    assert all(value>0 for value in signs)

    return {'sigma':case['sigma'],'delta':case['delta'],'phi':case['phi'],
            'points':{k:list(map(str,p)) for k,p in P.items()},'exact_distance_pairs':equal,
            'minimum_other_squared_distance_slack':str(min(slacks,key=lambda z:z.lower()).lower()),
            'minimum_ear_slack':str(min(ear_slacks,key=lambda z:z.lower()).lower()),
            'anchor_slacks':{k:str(value.lower()) for k,value in anchors.items()},
            'strict_m5_violation':str(violation.lower()),'center_translation':list(map(str,translation)),
            'F_unique_root_in_degrees':[str(left.lower()),str((phi+arb('0.3')).upper())],
            'F_min_positive_left_part_lower':str(min(signs,key=lambda z:z.lower()).lower()),
            'F_angle_derivative_per_radian_upper_on_neighborhood':str(derivative_upper.upper())}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    data=json.loads(args.input.read_bytes())
    results=[candidate(case) for case in data]
    args.output.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in row.items() if k!='points'} for row in results],indent=2))
