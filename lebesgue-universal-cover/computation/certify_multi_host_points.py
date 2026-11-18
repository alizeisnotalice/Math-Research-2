"""Arb replay of one finite diameter-one input against two hulls, with both hands."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import platform

import flint
from flint import arb, ctx

from certify_b9_padding_four_minimizer import axis, balance, features, SEED_DENOM
from certify_four_body_contact_stationary import BOX, add, det, dot, geometry, interval, normal, pi, q, r, regular

KSTAR_LOWER=F('.00315')
ESCAPE_LOWER=F('.00156')
LIMITS={'Kstar':KSTAR_LOWER,'R108':ESCAPE_LOWER}
SUPPORT_LIP={3:F('.29'),5:F('.163'),7:F('.115'),9:F('.089')}


class SupportHull:
    def __init__(self,pose,half=None):
        self.n=(3,5,7,9) if len(pose)==11 else (3,5,7)
        self.rot={3:arb(0),**{n:pose[k] for n,k in ((5,2),(7,5),(9,8)) if n in self.n}}
        offs={3:0,5:3,7:6,9:9}
        self.p={}
        self.hw={}
        h=[0]*len(pose) if half is None else half
        for n in self.n:
            angle=self.rot[n]
            k0=offs[n]
            v=normal(angle)
            t=pose[k0:k0+2]
            self.p[n]=[]
            self.hw[n]=(q(h[k0]),q(h[k0+1]),q(h[k0-1]) if n!=3 else arb(0))
            for k in range(n):
                z=regular(n,k)
                self.p[n].append(add(t,(v[0]*z[0]-v[1]*z[1],v[1]*z[0]+v[0]*z[1])))
        self.memo={}
        self.pads={n:q(SUPPORT_LIP[n]) for n in self.n}
        for n in self.n:
            assert (pi/(2*n)).tan()/2<self.pads[n]

    def direction(self,seed):
        if seed in self.memo:
            return self.memo[seed]
        w,a=axis(seed)
        u=tuple(map(q,w))
        h=r
        for n in self.n:
            h0=max((rad+dot(self.p[n][k],u)).upper() for k,rad,_ in features(n,a,self.rot[n]))
            ex,ey,ea=self.hw[n]
            upper=(h0+abs(u[0])*ex+abs(u[1])*ey+self.pads[n]*ea).upper()
            h=max(h,upper)
        row=(w,u,h)
        self.memo[seed]=row
        return row


def two_hulls(raw):
    data=json.loads(raw)
    K=geometry(*(interval(*s) for s in BOX))[0]
    R=[q(v) for v in data['zero_gain_five']['params']]
    return {'Kstar':SupportHull(K),'R108':SupportHull(R)}


def points(certificate):
    d=certificate['point_denominator']
    P=certificate['points']
    assert type(d) is int and d>0 and P[0]==[0,0]
    assert all(len(x)==2 and all(type(v) is int for v in x) for x in P)
    assert len({tuple(z) for z in P})==len(P)>=3
    D=max(sum((a-b)**2 for a,b in zip(x,y)) for x,y in combinations(P,2))
    assert D<d*d
    return [tuple(q(F(x,d)) for x in z) for z in P],D


def rotate(v,cs,sn):
    x,y=v
    return cs*x-sn*y,sn*x+cs*y


def lower(host,P,sg,dd,ii,facets):
    alpha=2*pi*q(F(2*ii+1,2**(dd+1)))
    w=pi/2**dd
    cs,sn=normal(alpha)
    normals=[host.direction(tuple(row[:2])) for row in facets]
    weights=balance([x[0] for x in normals])
    A=B=H=arb(0)
    for lam,(_,u,b),(_,_,index) in zip(weights,normals,facets,strict=True):
        assert type(index) is int and 0<=index<len(P)
        x,y=P[index]
        p=rotate((x,sg*y),cs,sn)
        A+=q(lam)*dot(p,u)
        B+=q(lam)*det(p,u)
        H+=q(lam)*b
    assert A>0 and 0<w<pi/2
    return (A*w.cos()-abs(B)*w.sin()-H).lower()


def verify_rows(host,P,sg,rows,target):
    partition=[]
    vals=[]
    for depth,index,facets in rows:
        assert type(depth) is int and 2<=depth<=45
        assert type(index) is int and 0<=index<2**depth
        assert len(facets)==3 and all(len(row)==3 for row in facets)
        interval=(F(index,2**depth),F(index+1,2**depth))
        partition.append(interval)
        lb=lower(host,P,sg,depth,index,facets)
        assert lb>q(target),(sg,depth,index,lb)
        vals.append(lb)
    partition.sort()
    assert partition[0][0]==0 and partition[-1][1]==1
    assert all(a[1]==b[0] for a,b in zip(partition,partition[1:]))
    return {'cells':len(rows),'area_independent_padding_lower':str(min(vals)),
            'depth':max(r[0] for r in rows),'allowed_padding':str(target)}


def replay(raw,certificate):
    assert certificate['competing_fits_input_sha256']==sha256(raw).hexdigest()
    assert certificate['padding_thresholds']=={k:str(v) for k,v in LIMITS.items()}
    hosts=two_hulls(raw)
    P,D=points(certificate)
    proofs={}
    for name,host in hosts.items():
        proofs[name]={}
        for s in (1,-1):
            sg=str(s)
            proofs[name][sg]=verify_rows(host,P,s,certificate['proofs'][name][sg],LIMITS[name])
    return {'python':platform.python_version(),'flint':flint.__version__,'precision':ctx.prec,
            'exact_pair_square_numerator':D,'max_pair_square_denominator':certificate['point_denominator']**2,
            'point_count':len(P),'tested_hands':proofs}


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('competing_fits',type=Path)
    parser.add_argument('certificate',type=Path)
    args=parser.parse_args()
    raw=args.competing_fits.read_bytes()
    cert=args.certificate.read_bytes()
    out=replay(raw,json.loads(cert))
    out['certificate_sha256']=sha256(cert).hexdigest()
    out['script_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
