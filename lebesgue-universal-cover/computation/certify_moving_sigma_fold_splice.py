"""Bound the changing ideal compatibility fold jointly in sigma and w."""
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys
from flint import arb_series

import certify_two_sigma_ideal_placement as b
from certify_low_sigma_fold_splice import ACTIVE, CYCLE

g=b.g
SIG=(F('1.54'),F('1.56'))
LEFT,RIGHT=F('52.65'),F('55.15')
W=(F('.42775'),F('.42957'))
P=F('.429')
TARGET=(F('1.5460'),F('1.5461'))


def D_of(sigma):
    e3=g.line_pair(b.directions[4],g.unit(7*g.PI/6+sigma))
    c=(-e3[0],e3[1])
    v=b.sub(b.z,c)
    vv=b.dot(v,v)
    ta=b.add(c,b.mul(1/vv,b.sub(v,b.mul((vv-1).sqrt(),b.rot(v)))))
    return b.sub(ta,b.z)


def harmonic(D,w):
    p=b.anchors(D,w)
    axes=[b.directions[i] for i in (0,2,4)]
    return tuple(sum(b.dot(f(u),p[j]) for u,(j,_) in zip(axes,ACTIVE,strict=True))
                 for f in (lambda u:u,b.rot))


def delta(D,w):
    a,c=harmonic(D,w)
    return a*a+c*c-g.HALF*g.HALF


def verify():
    lo,hi=SIG
    bot,top=W
    fl=delta(D_of(b.degree(hi)),g.exact(bot))
    fh=delta(D_of(b.degree(lo)),g.exact(top))
    assert fl>0 and fh<0
    target_lo,target_hi=(delta(D_of(b.degree(s)),g.exact(P)) for s in TARGET)
    assert target_lo>0 and target_hi<0,(target_lo,target_hi)
    counts=(64,256)
    cones=[]
    ds=[]
    dw=[]
    speeds=[]
    for i in range(counts[0]):
        s0,s1=[lo+(hi-lo)*F(k,counts[0]) for k in (i,i+1)]
        s=b.degree(s0).union(b.degree(s1))
        d_center=D_of(arb_series([b.degree((s0+s1)/2),g.PI/180],2))
        d_full=D_of(arb_series([s,g.PI/180,0],3))
        radius=b.ea((s0-s1)/2,(s1-s0)/2)
        speed=[m[1]+2*q[2]*radius for m,q in zip(d_center,d_full,strict=True)]
        D=[m[0]+v*radius for m,v in zip(d_center,speed,strict=True)]
        DD=[arb_series([x,v],2) for x,v in zip(D,speed,strict=True)]
        for j in range(counts[1]):
            w0,w1=[bot+(top-bot)*F(k,counts[1]) for k in (j,j+1)]
            w=b.ea(w0,w1)
            p=b.anchors(D,w)
            midangle=b.degree(LEFT).union(b.degree(RIGHT))
            for k,(idx,inv) in enumerate(ACTIVE):
                jv=CYCLE.index(idx)
                a,c=[b.sub(p[idx],p[CYCLE[(jv+r)%len(CYCLE)]]) for r in (-1,1)]
                if b.cross(a,c)<0:
                    a,c=c,a
                assert b.cross(a,c)>0
                u=g.unit(midangle+2*k*g.PI/3)
                if inv:
                    u=b.neg(u)
                cones += [b.cross(a,u).lower(),b.cross(u,c).lower()]
            a,c=harmonic(D,w)
            phi=g.arb.atan2(-c,-a)
            assert b.degree(LEFT)<phi<b.degree(RIGHT),(i,j,phi/b.g.PI*180,a,c)
            t1,t2=delta(DD,w)[1],delta(D,arb_series([w,1],2))[1]
            assert t1<0 and t2<0,(i,j,t1,t2)
            ds.append(t1);dw.append(t2);speeds.append(-t1/t2)
    assert min(cones)>0
    src=inspect.getsource(b.rig)
    target='scan(endcap,F(60)-endcap)'
    assert src.count(target)==1
    code=src.replace(target,'scan(endcap,F("52.65"));scan(F("55.15"),F(60)-endcap)').replace('assert b-a>min_cell','assert b-a>min_cell,(s0,s1,w0,w1,a,b,gap)')
    ns=b.rig.__globals__.copy()
    exec(compile(code,'rig_outside_moving_fold','exec'),ns)
    samples=(16,64)
    angular=0; outside=[]
    for i in range(samples[0]):
        s0,s1=[lo+(hi-lo)*F(k,samples[0]) for k in (i,i+1)]
        for j in range(samples[1]):
            w0,w1=[bot+(top-bot)*F(k,samples[1]) for k in (j,j+1)]
            nn,ee=ns['rig'](s0,s1,w0,w1,endcap=F('.3'),min_cell=F(1,1000000))
            angular+=nn
            outside.append(ee[-1])
    def enclose(v):
        return [str(min(a.lower() for a in v)),str(max(a.upper() for a in v))]
    return {'complete':True,'sigma_degrees':list(map(str,SIG)),'w':list(map(str,W)),
            'three_support_chamber_degrees':list(map(str,(LEFT,RIGHT))),
            'derivative_cells':list(counts),
            'Delta_w_on_all_cells':enclose(dw),'Delta_sigma_degree_on_all_cells':enclose(ds),
            'fold_w_sigma_degree_derivative':enclose(speeds),
            'sectors_min_cone_lower':str(min(cones)),
            'Delta_at_w_range_worst_corners':[str(fl),str(fh)],
            'pivot_fold_crossing_sigma_degrees':list(map(str,TARGET)),
            'Delta_at_pivot_crossing_sigma_ends':[str(target_lo),str(target_hi)],
            'angular_exclusion_cells':angular,'smallest_angular_exclusion_lower':str(min(outside)),
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':
    Path(sys.argv[1]).write_text(json.dumps(verify(),indent=2)+'\n')
