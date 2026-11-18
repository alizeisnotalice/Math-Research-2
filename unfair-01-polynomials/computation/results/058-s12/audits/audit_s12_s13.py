"""Independent exact S12 endpoint, field and explicit root-witness audit."""
from __future__ import annotations

from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import count
import json
from math import isqrt
from pathlib import Path
import sys
import time

import sympy as sp


def verify_formulas():
    x = sp.Symbol('x')
    a, aa, b, bb, c, d, e, f, g = sp.symbols('a abar b bbar c d e f g')
    m, u, r, s, omega = sp.symbols('m u r s omega', nonzero=True)
    gs = sp.symbols('G1:9')
    qs = sp.symbols('Q1:9')
    G=[1,*gs]
    Q=[1,*qs]
    P=[1,a,b,c,d,e,f,g,1]
    A5=sum(P[i]*G[5-i] for i in range(6))
    Pbar=[1,aa,bb,c,d,e,f,g,1]
    B5=sum(Pbar[i]*Q[5-i] for i in range(6))
    rho, sigma, t = sp.symbols('rho sigma t')
    subs={G[4]:rho-d,Q[4]:sigma-d,G[5]:-s*(rho-d)/r,Q[5]:-s*(sigma-d)/r}
    derivative=sp.simplify(sp.diff((A5-B5).subs(subs),d))
    assert sp.expand(derivative-(G[1]-Q[1]-a+aa))==0
    a1=a+G[1]
    identity=sp.expand(derivative.subs(Q[1],1-a-aa-G[1])-(2*a1-1-2*(a-aa)))
    assert identity==0
    aparts={a:m+r*omega,b:u+s*omega}
    for degree, sign in ((6,r*G[5]+s*G[4]),(7,r*G[6]+s*G[5])):
        value=sum(P[i]*G[degree-i] for i in range(degree+1))
        assert sp.diff(value.subs(aparts),omega)==sign

    # Symbolic distribution agrees for constant and short G and generic long G.
    for h in range(0,10):
        names=sp.symbols(f'j0:{h}')
        vals=([1]+list(names[1:]) if h else [])+[1]
        def get(k):
            return vals[k] if 0<=k<=h else 0
        def value(k,mother):
            return sum(coef*get(k-j) for j,coef in enumerate(mother))
        assert sp.expand(value(h+7,P)-(g+get(h-1)))==0
        assert sp.expand(value(h+6,P)-(f+g*get(h-1)+get(h-2)))==0
        assert sp.diff(value(h+2,P).subs(aparts),omega)==s
        assert sp.diff(value(h+1,P).subs(aparts),omega)==r+s*get(h-1)
        s13={a:m+r*omega,c:u+s*omega}
        assert sp.diff(value(h+3,P).subs(s13),omega)==s
        assert sp.diff(value(h+2,P).subs(s13),omega)==s*get(h-1)
        assert sp.diff(value(h+1,P).subs(s13),omega)==r+s*get(h-2)
    for degree, sign in ((5,r*G[4]+s*G[2]),(6,r*G[5]+s*G[3]),(7,r*G[6]+s*G[4])):
        val=sum(P[i]*G[degree-i] for i in range(degree+1))
        assert sp.diff(val.subs({a:m+r*omega,c:u+s*omega}),omega)==sign

    subs_a={a:sp.sqrt(3),b:sp.sqrt(3),c:5,d:5,e:6,g:1}
    minus2=sp.expand(sum(coef*(-2)**i for i,coef in enumerate(P)).subs(subs_a))
    assert minus2==64*f-23+2*sp.sqrt(3)
    return {'d_slope':str(derivative),'P_minus2':str(minus2),'endpoint_degrees':list(range(10))}


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def radical_sign(u,v,D):
    """Sign of u+v sqrt(D), by integer comparisons."""
    if not v:return (u>0)-(u<0)
    if not u or (u>0)==(v>0):return (v>0)-(v<0)
    q=u*u-D*v*v
    assert q
    return ((q>0)-(q<0))*((u>0)-(u<0))


def name_box(D, j):
    t=int(D%4==1)
    scale=1 if t else 2
    limit=(20 if j==1 else 83)
    mmax=(14 if j==1 else 60)
    kmax=(8 if j==1 else 37)
    assert scale*scale*D*(kmax+1)**2>=limit*limit
    rr=range(1,kmax+1) if j==1 else range(-kmax,kmax+1)
    result=set()
    for r in rr:
        if not r:continue
        assert 2*mmax+2-t*kmax>=limit
        for m in range(-mmax,mmax+1):
            T=2*m+t*r
            if radical_sign(limit-abs(T),-scale*abs(r),D)>0:
                result.add(((T,scale*r),r))
    return result


def verify_fields(source):
    sys.dont_write_bytecode=True
    spec=spec_from_file_location('frozen_s12_source',source)
    original=module_from_spec(spec)
    spec.loader.exec_module(original)
    fields=0
    raw=0
    extrema={}
    def record(key,value):
        old=extrema.setdefault(key,[value,value])
        old[:]=min(old[0],value),max(old[1],value)
    for D in range(2,400):
        if any(e!=1 for e in sp.factorint(D).values()):continue
        a=name_box(D,1)
        assert a==set(original.candidates(D,False)),('a',D)
        if not a:continue
        b=name_box(D,2)
        assert b==set(original.candidates(D,True)),('b',D)
        fields+=1
        raw+=len(a)*len(b)
        t=int(D%4==1)
        for label,rows in (('a',a),('b',b)):
            for (Ta,_),r in rows:
                assert (Ta-t*r)%2==0
                record(label+'_integer',(Ta-t*r)//2)
                record(label+'_omega',r)
    assert (fields,raw)==(112,1239790)
    return {'fields':fields,'raw_pairs':raw,'extrema':extrema}


def verify_rows(path, extra):
    rows=json.loads(extra.read_bytes())
    assert rows['reference_sha256']==sha(path)
    expected={tuple(row): (code,witness) for row,code,witness in rows['rows']}
    points=((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6))
    stats=Counter()
    for line in path.read_bytes().splitlines():
        row,code=json.loads(line)
        D,T,r,U,s,c,d,e,f,g=row
        assert r==s*g
        t=int(D%4==1)
        scale=1 if t else 2
        coeff=[(2,0),(T,scale*r),(U,scale*s),*((2*i,0) for i in (c,d,e,f,g,1))]
        def val(p,q):
            return tuple(sum(p**j*q**(8-j)*pair[i] for j,pair in enumerate(coeff)) for i in range(2))
        assert code==expected[tuple(row)][0]
        for p,q in points:
            u,v=val(p,q)
            for side, vv in enumerate((v,-v)):
                if radical_sign(u,vv,D)<=0:
                    break
            else:
                continue
            stats[f'positive_{p}/{q}']+=1
            witness=expected[tuple(row)][1]
            assert code=='root_positive_reject' and witness['point']==[p,q] and witness['side']==side
            assert witness['value_basis']==[(u-t*(vv//scale))//2,vv//scale]
            break
        else:
            assert code=='root_outer_reject'
            u,v=val(-2,1)
            assert radical_sign(u,v,D)<0
            stats['outer_negative_p_minus2']+=1
    assert stats=={'positive_1/1':13,'positive_2/1':2,'outer_negative_p_minus2':2}
    return {'rejections':dict(stats),'special_branch_mothers':sum(stats.values())}


def main():
    if not __debug__:raise SystemExit('requires assertion checking')
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    began=time.monotonic()
    path=args.source
    sc=path/'s12_exact_entry.py'
    primary=path/'full-entry.rows.jsonl'
    secondary=path/'independent-entry.json'
    result={'formulas':verify_formulas(),'basis':verify_fields(sc),'root_witnesses':verify_rows(primary,secondary),
            'python':sys.version,'sympy':sp.__version__,'seconds':time.monotonic()-began,
            'sources_sha256':{p.name:sha(p) for p in (sc,primary,secondary)},'status':'verified'}
    with args.out.open('x') as f:
        f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':
    main()
