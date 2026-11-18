"""Evaluate frozen balanced rows at all exact active edge-normal events."""
from argparse import ArgumentParser
from fractions import Fraction as F
from functools import cmp_to_key
from math import gcd
from pathlib import Path
import json
import sys

from flint import arb, ctx


def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]


def det(a,b):
    return a[0]*b[1]-a[1]*b[0]


def orientation_sort(a,b):
    def half(p):
        return 0 if p[1]>0 or (p[1]==0 and p[0]>0) else 1
    ha,hb=half(a),half(b)
    if ha!=hb:
        return ha-hb
    cross=det(a,b)
    return (cross<0)-(cross>0)


class Events:
    def __init__(self,points,den,hand,pi):
        from certify_nine import hull
        assert hand in (-1,1)
        self.points=hull([(x,hand*y) for x,y in points])
        assert len(self.points)>=3
        self.edges=[(b[1]-a[1],a[0]-b[0])
                    for a,b in zip(self.points,self.points[1:]+self.points[:1],strict=True)]
        self.den=den
        self.twopi=2*pi

    def endpoint(self,t,weights,n,den):
        c,s=self.cos_sin(t)
        total=arb(0)
        for w,(a,b) in zip(weights,n,strict=True):
            u,v=a*c+b*s,-a*s+b*c
            vals=[u*x+v*y for x,y in self.points]
            h=vals[0]
            for val in vals[1:]:
                h=h.max(val)
            total+=ar(w)*h
        return total/(den*self.den)

    def cos_sin(self,t):
        if t==0 or t==1:
            return arb(1),arb(0)
        z=self.twopi*ar(t)
        return z.cos(),z.sin()

    def event_value(self,c,s,weights,n,den):
        support=F(0)
        for w,(a,b) in zip(weights,n,strict=True):
            u,v=a*c+b*s,-a*s+b*c
            support+=w*max(u*x+v*y for x,y in self.points)
        assert support>=0
        return ar(support/(den*self.den))/arb(c*c+s*s).sqrt()

    def minimum(self,l,r,n,w,den):
        assert 0<=l<r<=1 and r-l<=F(1,2) and all(q>0 for q in w)
        assert all(sum(q*x[i] for q,x in zip(w,n,strict=True))==0 for i in (0,1))
        u,v=self.cos_sin(l),self.cos_sin(r)
        events=set()
        for m in self.edges:
            for a in n:
                c,s=dot(m,a),det(m,a)
                assert c or s
                if det(u,(c,s))<0 or det((c,s),v)<0:
                    continue
                g=gcd(abs(c),abs(s))
                events.add((c//g,s//g))
        events=sorted(events,key=cmp_to_key(orientation_sort))
        e=[self.endpoint(l,w,n,den),self.endpoint(r,w,n,den)]
        values=e+[self.event_value(c,s,w,n,den) for c,s in events]
        best=values[0]
        for val in values[1:]:
            best=best.min(val)
        return best,dict(events=[list(z) for z in events],endpoint_support=[str(t) for t in e],
                         minimum_support=str(best), event_support=[str(z) for z in values[2:]])


def read(p):
    return json.loads(p.read_bytes())


def cap_value(ev,caps,row):
    i,k,j,power,_=row
    l,r=F(i,2**k),F(i+1,2**k)
    needles,w,_=caps.triples[j]
    minimum,record=ev.minimum(l,r,[caps.host.normals[n] for n in needles],w,caps.host.den)
    bound=caps.bound(power,j)
    result=minimum-bound
    return result,record


def run(comp,one=False):
    sys.path.insert(0,str(comp))
    from certify_mixed_support_budget import Host, cloud
    from certify_high_six_host_neighborhood import Caps
    ctx.prec=160
    pi=arb.pi()
    data=read(comp/'high-six-host-neighborhood-angles.json')
    raw=read(comp/'multi-host-points.json')
    q,_=cloud(raw['points'],raw['point_denominator'])
    es={str(h):Events(raw['points'],raw['point_denominator'],h,pi) for h in (1,-1)}
    poses={
        'new_P30':read(comp/'free-p30-reoptimized-pose.json')[:11],
        'old_P30':read(comp/'free-p30-placement.json')[:11],
        'escape108':read(comp/'competing-fits-input.json')['zero_gain_five']['params'],
    }
    hosts={j:Host([F(a) for a in p],pi) for j,p in poses.items() if not one or j=='new_P30'}
    result={}
    for target,radii in data['config']['cap_flight']:
        for (host,d) in zip(poses,radii,strict=True):
            tag=f'{host}-{target}'
            if one and tag!='new_P30-.83755':
                continue
            c=Caps(hosts[host],q,pi,F(target),F(d))
            result[tag]={}
            for h,ev in es.items():
                if one and h!='1':
                    continue
                old=arb(1)
                new=arb(1)
                details=[]
                before=F(0)
                for j,row in enumerate(data['rows'][tag][h]):
                    i,k,index,n,pts=row
                    l,r=F(i,2**k),F(i+1,2**k)
                    assert l==before<r<=1
                    before=r
                    a=c.evaluate(l,r,int(h),index,n,pts)
                    b,record=cap_value(ev,c,row)
                    assert a>0 and b>0 and not b<a, (a,b,record)
                    old=old.min(a)
                    new=new.min(b)
                    record.update(row_index=j, angles=[str(l),str(r)],orig=str(a),exact_margin=str(b))
                    details.append(record)
                    if one:
                        break
                if not one:
                    assert before==1
                result[tag][h]=dict(cells=len(details),old_min=str(old),new_min=str(new),rows=details)
    return result


if __name__=='__main__':
    args=ArgumentParser(description=__doc__)
    args.add_argument('comp',type=Path)
    args.add_argument('--one',action='store_true')
    a=args.parse_args()
    print(json.dumps(run(a.comp,a.one),separators=(',',':')))
