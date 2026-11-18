"""Exact endpoint and algebraic-integrality certificates for C7=0."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import sympy as sp

t,u,w=sp.symbols('t u w')
H=t**6-4*t**5+6*t**4-7*t**3+4*t**2-2*t+1
m=1-t+t*t
n=1-t*m
p=1-(1-t)*n
v=t*t*m
q=1-n*v
x=1-p*v-(1-t)*q
y=1-p*v-t*q


def reduce(f):
    return sp.expand(sp.rem(sp.expand(f),H,t))


def interval(expr, lo, hi):
    out=(Q(0),Q(0))
    for a in sp.Poly(expr,t).all_coeffs():
        vals=tuple(z*s+Q(a) for z in out for s in (lo,hi))
        out=min(vals),max(vals)
    return out


def prefix():
    A=[1,t,0,n,p,0,0,0]+[0]*8
    B=[1,1-t,m,0,0,v,0,0]+[0]*8
    return A,B


def row(A,B,k):
    return sum(A[i]*B[k-i] for i in range(k+1))


def force(A,B,k,side):
    target=A if side=='A' else B
    assert target[k]==0
    target[k]=reduce(1-row(A,B,k))
    return target[k]


def cases():
    result={}
    # Degree 11. Locks A9, B8 and A10, B10.
    A,B=prefix()
    A[8]=q
    assert reduce(force(A,B,9,'B')-x)==0
    result['degree11_C10']=reduce(row(A,B,10))
    # Degree 12. Locks B8, B9, A10 and A11, B11.
    A,B=prefix()
    A[8]=q
    assert reduce(force(A,B,9,'A')-x)==0
    force(A,B,10,'B')
    result['degree12_C11']=reduce(row(A,B,11))
    # Degree 13. Locks A8, B9, B10, A11 and A12, B12.
    A,B=prefix()
    B[8]=q
    assert reduce(force(A,B,9,'A')-y)==0
    a10=force(A,B,10,'A')
    result['degree13_C12_lower']=reduce(p*q+m*a10)
    # Degree 14. Locks A9, B10, B11, A12 and A13, B13.
    A,B=prefix()
    A[8],B[8]=u,q-u
    force(A,B,9,'B')
    force(A,B,10,'A')
    force(A,B,11,'A')
    force(A,B,12,'B')
    expr=reduce(row(A,B,13)-1)
    assert sp.Poly(expr,u).degree()==1
    beta14=expr.coeff(u)
    gamma14=-expr.subs(u,0)
    assert sp.gcd(beta14,H)==1
    rbeta=int(sp.resultant(H,beta14,t))
    rgamma=int(sp.resultant(H,gamma14,t))
    norm=Q(rgamma,rbeta)
    assert norm==Q(8060939,208477) and norm.denominator>1
    # Degree 15. Locks A10, B11, B12, A13 and A14, B14.
    A,B=prefix()
    A[8],B[8]=u,q-u
    A[9]=w
    force(A,B,9,'B')
    B10=force(A,B,10,'B')
    A11=force(A,B,11,'A')
    force(A,B,12,'A')
    force(A,B,13,'B')
    eq15=reduce(row(A,B,14)-1)
    assert sp.Poly(eq15,u,w).total_degree()==1
    alpha15,beta15=eq15.coeff(u),eq15.coeff(w)
    gamma15=-eq15.subs({u:0,w:0})
    assert beta14==beta15
    assert reduce(B10.subs({u:q,w:0})-(1-m*q-t*x))==0
    at_zero=reduce(A11.subs({u:q,w:x}, simultaneous=True))
    assert at_zero==reduce(1-t*(1-m*q-(1-t)*x)-m*x)
    vals=[int(f.subs(t,-6))%163 for f in (H,alpha15,beta15,gamma15)]
    assert vals==[0,0,0,132]
    result['degree15_A11_if_C12_zero']=at_zero
    return result,{
        'degree14_beta':str(beta14), 'degree14_gamma':str(gamma14),
        'degree14_norm_beta':rbeta,'degree14_norm_gamma':rgamma,'degree14_norm_A8':str(norm),
        'degree15_alpha':str(alpha15),'degree15_beta':str(beta15),'degree15_gamma':str(gamma15),
        'degree15_modulus':163,'degree15_t_residue':-6,'degree15_residues_H_alpha_beta_gamma':vals,
    }


def main():
    if not __debug__:raise SystemExit('run without -O')
    options=ArgumentParser(description=__doc__)
    options.add_argument('--out',required=True,type=Path)
    o=options.parse_args()
    start=time.monotonic()
    sturm=sp.sturm(H,t)
    signs=[[int(sp.sign(g.subs(t,z))) for g in sturm] for z in (0,1)]
    def variation(signs):
        active=[s for s in signs if s]
        return sum(a!=b for a,b in zip(active,active[1:]))
    variations=list(map(variation,signs))
    assert variations[0]-variations[1]==1
    lo,hi=Q(71942965,10**8),Q(71942966,10**8)
    assert H.subs(t,lo)>0>H.subs(t,hi)
    polys,data=cases()
    polys.update(q=q,x=x,y=y,minimum_10_for_14_and_15=1-m*q-t*x)
    boxes={}
    wanted={
        'q':(Q('0.82410'),Q('0.82411')),'x':(Q('0.40502'),Q('0.40504')),
        'y':(Q('0.04335'),Q('0.04337')),'minimum_10_for_14_and_15':(Q('0.05085'),Q('0.05086')),
        'degree11_C10':(Q('0.9491'),Q('0.9492')),
        'degree12_C11':(Q('0.4877'),Q('0.4878')),
        'degree13_C12_lower':(Q('1.5140'),Q('1.5142')),
        'degree15_A11_if_C12_zero':(Q('0.5122'),Q('0.5123')),
    }
    for key,expr in polys.items():
        red=reduce(expr)
        bounds=interval(red,lo,hi)
        assert wanted[key][0]<bounds[0]<=bounds[1]<wanted[key][1],(key,tuple(map(float,bounds)))
        boxes[key]={'poly':str(red),'lo':str(bounds[0]),'hi':str(bounds[1]),
                    'display_lo':str(wanted[key][0]),'display_hi':str(wanted[key][1])}
    data.update(status='verified', sturm_endpoint_signs=signs, sturm_variations=variations,
                t_isolating_interval=[str(lo),str(hi)], verified_boxes=boxes,python=sys.version,
                sympy=sp.__version__,time_seconds=time.monotonic()-start,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    with o.out.open('x') as f:f.write(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:data[k] for k in ('status','degree14_norm_A8','degree15_residues_H_alpha_beta_gamma','time_seconds')},sort_keys=True))


if __name__=='__main__':main()
