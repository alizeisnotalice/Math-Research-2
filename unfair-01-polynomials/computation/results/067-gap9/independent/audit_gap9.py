#!/usr/bin/env python3
"""Independently audit ninth-row branches with lex elimination and rational Sturm."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp
from prefix_gap import T, equations, primitive, shapes, strip


def trim(p):
    while p and not p[-1]:p.pop()
    return p


def rem(a,b):
    a=a[:]
    while len(a)>=len(b):
        d=len(a)-len(b)
        q=a[-1]/b[-1]
        for j,x in enumerate(b):a[d+j]-=q*x
        trim(a)
    return a


def at(p,x):
    out=Q(0)
    for a in reversed(p):out=out*x+a
    return out


def coeffs(expr):
    return [Q(x) for x in sp.Poly(expr,T,domain=sp.QQ).all_coeffs()[::-1]]


def sturm(poly):
    c=coeffs(poly)
    out=[c,[i*c[i] for i in range(1,len(c))]]
    while out[-1]:
        r=rem(*out[-2:])
        if not r:break
        out.append([-x for x in r])
    return out


def variation(chain,x):
    v=[(u>0)-(u<0) for p in chain if (u:=at(p,x))]
    return sum(a!=b for a,b in zip(v,v[1:]))


def number_open(p):
    chain=sturm(p)
    return variation(chain,Q(0))-variation(chain,Q(1))-int(p.eval(1)==0)


def isolate(p):
    chain=sturm(p)
    lo,hi=Q(0),Q(1)
    assert variation(chain,lo)-variation(chain,hi)==1
    for _ in range(64):
        mid=(lo+hi)/2
        assert at(chain[0],mid)
        if variation(chain,lo)>variation(chain,mid):
            hi=mid
        else:
            lo=mid
    assert 0<lo<hi<1
    return lo,hi


def interval(expr,ends):
    lo,hi=ends
    a=b=Q(0)
    for c in reversed(coeffs(expr)):
        values=(lo*a,lo*b,hi*a,hi*b)
        a,b=min(values)+c,max(values)+c
    return a,b


def check_sign(shape,p,raw,A,B,symbols):
    x,t=symbols
    assert t==T and len(raw)==2
    f,g=map(lambda r:sp.Poly(r,x),raw)
    relation=sp.Poly(sp.expand(f.as_expr()*g.LC()-g.as_expr()*f.LC()),x)
    assert relation.degree()==1
    slope,offset=relation.all_coeffs()
    if shape=="AXB":
        numerator=-offset
    else:
        assert shape=="BXA"
        numer=sp.Poly(B[3],x)
        assert numer.nth(1)==-1
        numerator=sp.expand(numer.nth(0)*slope+offset)
    lo,hi=isolate(p)
    U=interval(slope,(lo,hi))
    X=interval(numerator,(lo,hi))
    assert U[0]*U[1]>0
    values=[n/d for n in X for d in U]
    bound=(min(values),max(values))
    assert bound[1]<Q(-2,5)
    return dict(shape=shape,isolating_interval=list(map(str,(lo,hi))),
                coefficient="A3" if shape=="AXB" else "B3",
                bound=list(map(str,bound)),denominator=[str(v) for v in U])


def main():
    if not __debug__:raise SystemExit("assertions required")
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    raw=args.reference.read_bytes()
    original=json.loads(raw)
    assert [r["shape"] for r in original["ledger"]]==list(shapes())
    proofs=[]
    factors_by_shape={}
    for row in original["ledger"]:
        shape=row["shape"]
        if shape[:2]=="XX":
            assert row["reason"]=="eighth_row_at_most_one_product"
            continue
        if shape[0]==shape[2]=="X":
            assert row["reason"]=="eighth_row_quadratic_upper_bound"
            continue
        symbols,raw_eq,A,B=equations(9,shape)
        primitive_eq=[strip(f,symbols) for f in raw_eq]
        G=sp.groebner(primitive_eq,*symbols,domain=sp.QQ,order="lex")
        eliminants=[f for f in G if f.free_symbols<={T}]
        assert len(eliminants)==1 and eliminants[0]
        independent=primitive(eliminants[0])
        product=sp.Poly(1,T,domain=sp.QQ)
        for item in row["factors"]:
            p0=sp.Poly.from_list(item["coefficients"],T,domain=sp.QQ)
            assert p0.is_irreducible
            assert number_open(p0)==item["open_roots"]
            product*=p0**item["exponent"]
        # Resultants are ideal consequences. Endpoint factors may be normalized away.
        check=primitive(strip(independent.as_expr(),(T,)))
        assert product.rem(check).is_zero,(shape,product,check)
        case=dict(shape=shape,independent_degree=independent.degree(),signs=[])
        for f,_ in sp.factor_list(check.as_expr())[1]:
            p=primitive(f)
            assert p.is_irreducible
            if p.LC()!=1 or number_open(p)==0:continue
            assert p.degree()==16
            hit=check_sign(shape,p,raw_eq,A,B,symbols)
            factors_by_shape[shape]=p
            case["signs"].append(hit)
        proofs.append(case)
    assert set(factors_by_shape)=={"AXB","BXA"}
    p,q=factors_by_shape["AXB"],factors_by_shape["BXA"]
    assert sp.expand(p.as_expr().subs(T,1-T))==q.as_expr()
    result=dict(status="verified",checked_symbolic=len(proofs),checked_analytic=27-len(proofs),
                reference_sha256=sha256(raw).hexdigest(),ledger=proofs)
    with args.out.open("x") as f:f.write(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("status","checked_symbolic","checked_analytic")}))


if __name__=="__main__":main()
