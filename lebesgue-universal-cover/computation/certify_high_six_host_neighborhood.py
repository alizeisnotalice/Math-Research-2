"""Joint long-edge caps and inherited mixed-area rows above the six-body upper example."""
from argparse import ArgumentParser
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
from pathlib import Path
import json
import sys

from flint import arb, ctx

from certify_mixed_support_budget import (
    CENTERED_P30_POSES, Host, Problem, PI_UP, ar, cloud, det, dot, read, rotate, sha,
)

LOW = F('.83740')
R = F('.703')
CAP_FLIGHT = (
    ('.83753', ('.00022', '.00025', '.00028')),
    ('.83755', ('.00013', '.00017', '.00018')),
)
MIXED_TARGETS = (F('.83760'), F('.83766'))
KEYS = ('new_P30','old_P30','escape108')
MODES = (1, 2, 4, 8, 12, 24)


def edge_lengths(host):
    return [arb(sq).sqrt()/host.den for sq, _, _ in host.edges]


def perimeter_upper(host):
    poly = host.material
    length = 0
    for a,b in zip(poly,poly[1:]+poly[:1],strict=True):
        qq = (a[0]-b[0])**2 + (a[1]-b[1])**2
        r = isqrt(qq)
        length += r+(r*r<qq)
    return F(length,host.den)


def groups_for(normals):
    acute = [(a,b) for a,b in combinations(range(3),2) if dot(normals[a],normals[b])>0]
    assert len(acute)<=1
    if not acute:
        return ((0,),(1,),(2,))
    group=acute[0]
    return group,tuple(j for j in range(3) if j not in group)


class Caps:
    def __init__(self,host,pts,pi,target,radius):
        self.pi,self.host,self.pts=pi,host,pts
        self.delta=target-LOW
        self.d=radius
        self.s=1-2*radius
        self.lens=edge_lengths(host)
        self.base=host.offsets('mixed')
        self.depth=[ar(radius+host.eta)+(1-ar(self.s))*b/ell
                    for b,ell in zip(self.base,self.lens,strict=True)]
        self.edge_limit=[b+ell*ar(radius+host.eta) for b,ell in zip(self.base,self.lens,strict=True)]
        cap_budget=(ar(self.delta)*(4*ar(R)**2-1).sqrt()).sqrt()
        enclosing=ar(F(3,5)+radius+host.eta)+cap_budget
        assert enclosing<ar(F('.613'))
        assert enclosing<arb(2).sqrt()*ar(F(1,2)-radius)
        self.enclosing=str(enclosing)
        admissible=[i for i,(sq,_,_) in enumerate(host.edges) if sq*1600>host.den**2]
        self.triples=[]
        for ix in combinations(admissible,3):
            normals=[host.normals[i] for i in ix]
            a,b,c=normals
            w=(det(b,c),det(c,a),det(a,b))
            if all(x<0 for x in w):
                w=tuple(-x for x in w)
            if not all(x>0 for x in w):
                continue
            mass=sum(w)
            coefficients=tuple(F(x,mass) for x in w)
            assert all(sum(w*n[j] for w,n in zip(coefficients,normals,strict=True))==0 for j in (0,1))
            self.triples.append((ix, coefficients, groups_for(normals)))
        self.cache={}

    def bound(self,power,index):
        k=(power,index)
        if k in self.cache:
            return self.cache[k]
        assert power in MODES
        ids,w,groups=self.triples[index]
        masses=[sum(w[i] for i in group)**power for group in groups]
        linear=arb(0)
        slope=arb(0)
        for group,mass in zip(groups,masses,strict=True):
            a=ar(self.delta*mass/sum(masses))
            g=arb(0)
            for i in group:
                j=ids[i]
                b=ar(self.s)*self.lens[j]*self.depth[j]
                root=(a*a+2*b*a).sqrt()
                v=ar(w[i]/self.s)*(a+root)
                prime=ar(w[i]/self.s)*(1+(a+b)/root)
                linear+=v-a*prime
                g+=prime
            slope=slope.max(g)
        value=linear+slope*ar(self.delta)
        value+=sum((ar(k)*self.edge_limit[j] for j,k in zip(ids,w,strict=True)),arb(0))
        self.cache[k]=value
        return value

    def evaluate(self,left,right,hand,index,power,points):
        angle=ar(left+right)*self.pi
        assert 0<right-left<=F(1,2)
        width=ar(right-left)*self.pi
        c,s=angle.cos(),angle.sin()
        ids,w,_=self.triples[index]
        assert len(points)==3
        A=B=arb(0)
        for j,r,k in zip(ids,w,points,strict=True):
            assert 0<=k<len(self.pts)
            x,y=map(ar,self.pts[k])
            p=rotate((x,y*hand),c,s)
            n=self.host.nums[j]
            A+=ar(r)*dot(p,n)
            B+=ar(r)*det(p,n)
        assert A>0
        return A*width.cos()-abs(B)*width.sin()-self.bound(power,index)

    def make(self,hand):
        import numpy as np

        n=np.array(self.host.normals,dtype=float)/self.host.den
        cloud=np.array(self.pts,dtype=float)*[1,hand]
        weights=np.zeros((len(self.triples),len(n)))
        shifts=np.zeros(len(self.triples))
        options=[]
        for j,(idx,w,_) in enumerate(self.triples):
            weights[j,list(idx)]=list(map(float,w))
            val,power=min((float(self.bound(power,j).mid()),power) for power in MODES)
            options.append(power)
            shifts[j]=val

        def cell(i,k):
            l,r=F(i,2**k),F(i+1,2**k)
            ang=float(l+r)*np.pi
            c,s=np.cos(ang),np.sin(ang)
            support=n@(cloud@np.array([[c,s],[-s,c]])).T
            best=support.argmax(axis=1)
            gaps=weights@support.max(axis=1)-shifts
            for j in np.argsort(gaps)[-12:][::-1]:
                idx=int(j)
                pts=[int(best[q]) for q in self.triples[idx][0]]
                p=options[idx]
                if self.evaluate(l,r,hand,idx,p,pts)>0:
                    yield [i,k,idx,p,pts]
                    return
            assert k<17,(hand,i,k,float(gaps.max()))
            yield from cell(2*i,k+1)
            yield from cell(2*i+1,k+1)
        return [row for i in range(64) for row in cell(i,6)]

    def check(self,hand,rows):
        prev=F(0)
        smallest=arb(1)
        for i,k,j,p,pts in rows:
            assert 0<=j<len(self.triples)
            l,r=F(i,2**k),F(i+1,2**k)
            assert prev==l<r<=1
            value=self.evaluate(l,r,hand,j,p,pts)
            assert value>0
            smallest=smallest.min(value)
            prev=r
        assert prev==1
        return dict(cells=len(rows),minimum_margin=str(smallest))


def mixed_radius(host,pts,pi,old,key,target):
    """Reuse the original mixed-area rows with a directed Hausdorff budget."""
    tag=f'{key}-Q11-mixed'
    cells=old['proofs'][tag]
    problem=Problem(host,pts,'mixed',pi,F(0))
    perimeter=perimeter_upper(host)
    limit=arb(1)
    for sign in (1,-1):
        prev=F(0)
        for i,k,rows in cells[str(sign)]:
            l,r=F(i,2**k),F(i+1,2**k)
            assert l==prev<r<=1
            g=problem.evaluate(l,r,sign,rows)
            L=sum((ar(F(a,b))*(arb(host.edges[j][0]).sqrt()/host.den)
                   for j,_,a,b in rows),arb(0))
            coefficient=ar(perimeter/2+2*target)-L
            assert coefficient>0
            window=(ar(host.area-target)+g)/coefficient
            assert window>0
            limit=limit.min(window)
            prev=r
        assert prev==1
    radius=F(int(float(limit.mid())*10**7)-1,10**7)
    assert 0<ar(radius)<limit
    return dict(radius=str(radius),limit_enclosure=str(limit))


def run(a):
    ctx.prec=160
    sys.path.insert(0,str(a.computation))
    files={k:a.computation/(name+'.json') for k,name in {
        'new_P30':'free-p30-reoptimized-pose','old_P30':'free-p30-placement',
        'escape108':'competing-fits-input','Q11':'multi-host-points','P30':'finite-point-leaf-input',
    }.items()}
    inputs={k:sha(p) for k,p in files.items()}
    raw={k:read(p) for k,p in files.items()}
    old=read(a.mixed_angles)
    assert all(inputs[k]==old['config']['inputs'][k] for k in inputs)
    poses={k:raw[k][:11] for k in KEYS[:2]}
    poses[KEYS[2]]=raw[KEYS[2]]['zero_gain_five']['params']
    pi=arb.pi()
    assert pi<ar(PI_UP)
    def area_point(radius):
        x=ar(radius)
        return pi/4+(x*x-ar(F(1,4))).sqrt()/2-(1/(2*x)).acos()/4

    assert area_point(R)>ar(max(F(r[0]) for r in CAP_FLIGHT))
    mixed_disk_bound=F('.704')
    assert area_point(mixed_disk_bound)>ar(max(MIXED_TARGETS))
    hosts={k:Host(list(map(F,poses[k])),pi) for k in KEYS}
    assert all(h.hull_hash==old['config']['host_inner'][k] for k,h in hosts.items())
    q11,q11_d=cloud(raw['Q11']['points'],raw['Q11']['point_denominator'])
    p30,p30_d=cloud(raw['P30']['points'],raw['P30']['denominator'])
    config=dict(inputs=inputs,original_angles_sha256=sha(a.mixed_angles),
                cap_flight=CAP_FLIGHT,area_lower=str(LOW),disk_bound=str(R))
    # JSON stores the tuple-valued constants as lists.
    config=json.loads(json.dumps(config))
    if a.generate and not a.certificate.exists():
        recorded={}
    else:
        data=read(a.certificate)
        assert data['config']==config
        recorded=data['rows']
    result={}
    for t,radii in CAP_FLIGHT:
        for key,d in zip(KEYS,radii,strict=True):
            problem=Caps(hosts[key],q11,pi,F(t),F(d))
            tag=f'{key}-{t}'
            if a.generate and tag not in recorded:
                recorded[tag]={str(s):problem.make(s) for s in (1,-1)}
                a.certificate.write_text(json.dumps(dict(config=config,rows=recorded),separators=(',',':'))+'\n')
            summary={str(s):problem.check(s,recorded[tag][str(s)]) for s in (1,-1)}
            result[tag]=dict(radius=str(F(d)),results=summary,anchor_frame_radius=problem.enclosing)
            if a.generate:
                print(tag,d,summary,file=sys.stderr,flush=True)
    high={}
    for key,h in hosts.items():
        # Both curve-hull area and a true P30 fit are checked, so these anchors
        # actually lie in the six-body sublevel for the top transfer.
        area=h.area+perimeter_upper(h)*h.eta+PI_UP*h.eta**2
        assert area<F('.83766')
        fit=h.fits(p30,CENTERED_P30_POSES[key],pi)
        bounds={str(t):mixed_radius(h,q11,pi,old,key,t) for t in MIXED_TARGETS}
        high[key]=dict(witness_fit_clearance=fit,anchor_area_upper=str(area),mixed_bounds=bounds)
    return dict(config=config,cap=result,higher=high,diameter_squares=[q11_d,p30_d],
                certificate_sha256=sha(a.certificate),source_sha256=sha(Path(__file__)),mixed_disk_bound=str(mixed_disk_bound))


if __name__=='__main__':
    arg=ArgumentParser(description=__doc__)
    for name in ('computation','mixed_angles','certificate'):
        arg.add_argument(name,type=Path)
    arg.add_argument('--generate',action='store_true')
    print(json.dumps(run(arg.parse_args()),separators=(',',':')))
