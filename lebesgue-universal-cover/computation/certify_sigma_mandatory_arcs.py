"""Verify three compulsory disk intersections over a closed, moving P-host angle interval."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from audit_eh_continuous_placement import project
from certify_eh_general_relocation import rotate_angle


def degree(f):
    return g.exact(f)*g.PI/180


def anchors(label):
    n,m,_,b,_,_,s,t,_=g.setup()
    if label=='TI':
        axis=5
        r=g.line_pair(m[2],n[1],g.HALF,g.arb(0))
        f=g.unit(5*g.PI/6-g.SIGMA)
        a=g.scale(g.arb('.998'),tuple(map(g.arb,['.125327173657','.467705236099'])))
        pts={'S':s,'E':g.minus(s,n[0]),'R':r,'U':g.minus(r,f),'B':b,'Z':g.minus(b,n[1]),'T':t,'a':a}
        copies=('S','E','U','T','a')
        pairs=[('S','E'),('S','R'),('S','T'),('R','U'),('B','Z'),('B','T'),
               ('R','JU'),('R','JS'),('JS','JT'),('JS','JE'),('Z','JT')]
        phase,three,guards=1,('E','JS','B'),[(1,'JT')]
    elif label=='GK':
        axis=1
        c=g.line_pair(n[3],m[2])
        pts={'C':c,'E':g.plus(c,n[0]),'G':g.circle_line(c,n[5],1),'K':project(c,m[5]),
             'A':g.scale(g.arb('.56'),g.unit(g.PI/6))}
        copies=('C','E','G','K')
        pairs=[('C','E'),('C','G'),('C','K'),('G','JG'),('JC','JE'),('JC','JG'),('JC','JK')]
        phase,three,guards=0,('E','JG','JC'),[(0,'A'),(1,'JG')]
    else:
        axis=3
        c=g.line_pair(n[4],m[3])
        pts={'C':c,'Q':project(c,n[1]),'L':g.plus(c,m[0]),'Y':(g.arb(0),c[1]+(1-c[0]**2).sqrt()),
             'A':(g.HALF,g.arb('.038'))}
        copies=('C','Q','L','A')
        pairs=[('C','Q'),('C','L'),('C','Y'),('Y','JC'),('A','JA'),('JC','JQ'),('JC','JL')]
        phase,three,guards=1,('JA','JC','Q'),[(1,'Y'),(2,'JL')]
    mirror=lambda z:rotate_angle((z[0],-z[1]),axis*g.PI/3)
    pts.update({'J'+k:mirror(pts[k]) for k in copies})
    return pts,pairs,phase,three,guards


def certify_cell(label,left,right):
    g.SIGMA=degree(left).union(degree(right))
    pts,pairs,phase,three,guards=anchors(label)
    _,m,_,b,_,_,s,t,_=g.setup()
    skip={frozenset(x) for x in pairs}
    slacks=[(1-g.dot(g.minus(v,w),g.minus(v,w))).lower() for (i,v),(j,w) in combinations(pts.items(),2)
            if frozenset((i,j)) not in skip]
    assert min(slacks)>0
    lam=g.arb('.0000001')
    def axes(theta):
        return [g.unit(phase*g.PI+theta+2*k*g.PI/3) for k in range(3)]
    B=sum(g.dot(g.rot(u),pts[k]) for u,k in zip(axes(g.arb(0)),three,strict=True))
    arc1=degree(F(1))
    deriv=(1-lam)*(B*arc1.cos()-3*g.HALF*arc1.sin())-9*lam/100
    side=(3*g.SIGMA).cos()/100
    caps=min(((1-lam)*(g.dot(m[i],pts[k])-g.HALF)-lam/100).lower() for i,k in guards)
    assert min(deriv,side,caps)>0
    signs=[]
    def partition(a,b):
        theta=degree(a).union(degree(b))
        floor=sum(max(g.dot(u,v).lower() for v in pts.values()) for u in axes(theta))-3*g.HALF
        score=(1-lam)*floor-3*lam/100
        if score>0:
            signs.append(score.lower())
            return
        assert b-a>F(1,1024)
        mid=(a+b)/2
        partition(a,mid)
        partition(mid,b)
    partition(F(1),F(30))
    stats={'diameter_slack':min(slacks),'root_slope':deriv,'strict_caps':caps,
           'side_excess_per_weight':side,'interior_root_excess':min(signs)}
    def arc(v,w,key):
        assert g.cross(v,w)>0 and g.dot(v,w)>0
        ang=g.arb.atan2(g.cross(v,w),g.dot(v,w))
        stats[key]=ang.lower()
        stats[key+'_negative_upper']=-ang.upper()
        return ang
    if label=='TI':
        safe_radius=g.arb('.04')
        margin=min((1-safe_radius-g.norm(g.minus(t,z))).lower() for k,z in pts.items() if k not in ('S','B'))
        rad=g.minus(t,b)
        ang=arc(rad,g.neg(m[0]),'TI_radians')
        deriv2=2*g.dot(g.minus(b,s),g.rot(rotate_angle(rad,g.arb(0).union(ang))))
        assert margin>0 and ang<safe_radius and deriv2<0
        stats['lens_distance_slack']=margin
        stats['minus_arc_distance_derivative']=-deriv2.upper()
    elif label=='GK':
        arc(g.minus(pts['G'],pts['C']),g.minus(pts['K'],pts['C']),'GK_radians')
    else:
        c=pts['C']
        y=pts['Y']
        arc(g.minus(pts['Q'],c),g.minus(y,c),'QY_radians')
        jc=pts['JC']
        rt=g.arb(3).sqrt()
        z=g.circle_pair((rt/2-g.HALF,(rt/2-g.HALF-1)/rt),(-g.HALF/2,-rt/4),-1)
        v=g.minus(z,jc)
        vv=g.dot(v,v)
        ta=g.plus(jc,g.scale(1/vv,g.minus(v,g.scale((vv-1).sqrt(),g.rot(v)))))
        u0,u1,u2=[g.minus(j,jc) for j in (y,ta,pts['JQ'])]
        assert g.cross(u1,u2)>0 and g.dot(u0,u2)>0
        arc(u0,u1,'YTa_radians')
    return len(signs),{k:z.lower() for k,z in stats.items()}


def verify(left,right):
    assert F(1)<=left<right<=F(3)
    result={}
    for label in ('TI','GK','QY'):
        cells=[]
        def visit(a,b):
            try:
                count,stats=certify_cell(label,a,b)
            except AssertionError:
                assert b-a>F(1,20000),(label,a,b)
                mid=(a+b)/2
                visit(a,mid)
                visit(mid,b)
                return
            cells.append((a,b,count,stats))
        visit(left,right)
        result[label]={'sigma_parts':[[str(a),str(b),n] for a,b,n,_ in cells],
                       'rectangle_count':sum(n for _,_,n,_ in cells),'lower_bounds':{
                           k:str(min(stats[k] for _,_,_,stats in cells)) for k in cells[0][3]}}
    return {'sigma_degrees':[str(left),str(right)],'root_degrees':['0','60'],
            'precision_bits':g.ctx.prec,'cases':result}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('low',type=F)
    p.add_argument('high',type=F)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(verify(args.low,args.high),indent=2)+'\n')
