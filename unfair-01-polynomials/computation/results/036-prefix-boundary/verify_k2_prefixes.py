#!/usr/bin/env python3
"""Exact Fraction audit of the two K=2 q=21 prefix witnesses.

No search module is imported.  Each supplied F prefix is separately extended
by zero only for finite polynomial checks; this is deliberately not a claim
about a common full F or a full factorization.
"""
from __future__ import annotations

import hashlib, json, sys, time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INPUTS = (ROOT/'k2-prefix.json', ROOT/'k2-reverse-prefix.json')

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def conv(a,b,degree):
    return [sum((a[i]*b[d-i] for i in range(len(a)) if 0 <= d-i < len(b)), Fraction(0)) for d in range(degree+1)]
def compose_product(a,u,b,v):
    out=[Fraction(0)]*(u*(len(a)-1)+v*(len(b)-1)+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[u*i+v*j]+=x*y
    return out

def audit(path):
    data=json.loads(path.read_text()); p,q,s,t=data['p'],data['q'],data['s'],data['t']
    if (p,q,t,data['prefix_length']) != (2,21,20,42) or s not in (8,12): raise ValueError('unexpected witness parameters')
    a=Fraction(data['a']); F=[Fraction(x) for x in data['first']['F_prefix']]; H=[Fraction(x) for x in data['first']['H_prefix']]
    if a != Fraction(1,2) or len(F)!=42 or len(H)!=42 or F[0]!=1 or any(x not in (0,1) for x in H): raise ValueError('bad exact prefix encoding')
    G=[Fraction(0)]*42
    for d,c in ((0,1),(s,a),(t,1),(q,1),(q+s,a),(q+t,1)):
        if G[d]: raise ValueError('G support collision')
        G[d]=Fraction(c)
    fg=conv(F,G,82)
    if fg[:42] != H: raise ValueError('FG does not reproduce supplied H prefix')
    if any(x<0 for x in F): raise ValueError('negative F prefix')
    # The prescribed finite crossed rows.  The code intentionally checks only
    # these intervals, retaining all later coefficients as out of scope.
    A=compose_product(F,p,G,q); B=compose_product(F,q,G,p)
    if len(A)!=944 or len(B)!=944: raise AssertionError('degree arithmetic')
    if not all(0<=x<=1 for x in A[:84]): raise ValueError('A prefix outside [0,1]')
    if not all(0<=x<=1 for x in B[:882]): raise ValueError('B prefix outside [0,1]')
    # Hhat is precisely this length-42 H with no asserted tail.  Its complete
    # composition product is nevertheless a finite, exact polynomial.
    C=compose_product(H,p,H,q)
    if any(x not in (0,1) for x in C): raise ValueError('Hhat(x^2)Hhat(x^21) is not binary')
    first_tail=next((d for d in range(42,len(fg)) if fg[d] != 0),None)
    first_nonbinary=next((d for d in range(42,len(fg)) if fg[d] not in (0,1)),None)
    expected={8:(44,Fraction(5,8)),12:(48,Fraction(5,16))}[s]
    if first_nonbinary != expected[0] or fg[first_nonbinary] != expected[1]: raise ValueError('unexpected first truncated FG failure')
    return {'source':str(path),'source_sha256':sha(path),'s':s,'a':str(a),'G_support':[i for i,x in enumerate(G) if x],
            'F_nonnegative_through':41,'FG_matches_H_through':41,
            'A_checked_through':83,'B_checked_through':881,
            'crossed_prefix_ranges':'all coefficients in [0,1]',
            'Hhat_product_degree':len(C)-1,'Hhat_product_binary':True,
            'first_FG_tail_discrepancy_degree':first_tail,'first_FG_tail_discrepancy_value':str(fg[first_tail]),'first_full_product_nonbinary_degree':first_nonbinary,'first_full_product_nonbinary_value':str(fg[first_nonbinary])}

def main():
    started=time.monotonic(); results=[audit(x) for x in INPUTS]
    out=Path(__file__).with_name('k2-prefix-audit-receipt.json')
    rec={'audit':'independent-k2-prefix-fraction-v1','scope':'Two separate finite prefix witnesses only: neither supplies a full F tail, a full factorization, nor a counterexample.','results':results,
         'audit_source_sha256':sha(Path(__file__)),'interpreter':sys.version,'elapsed_seconds':time.monotonic()-started,'status':'complete'}
    out.write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n'); print(json.dumps({'results':results,'elapsed_seconds':rec['elapsed_seconds'],'status':'complete'},sort_keys=True))
if __name__=='__main__':main()
