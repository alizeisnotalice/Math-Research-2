#!/usr/bin/env python3
"""Exact root-by-root necessary-prefix audit for 3<=r<=8.

For each endpoint branch u=c_(2r), isolate every zero in 0<s<1 of its
required residual (F for u=1, F+s for u=0), then record a failed strict
condition t in (0,1) or a negative prefix coefficient.  All root counts use
SymPy QQ Sturm arithmetic; the report deliberately retains rational boxes.
"""
import json, os
from pathlib import Path
import sympy as S
s=S.symbols('s'); LO=S.Rational(0); HI=S.Rational(1)
out=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round011/generated-small-ledger.json')
rows=[]
for r in range(3,9):
    D=1-2*s+3*s*s; U=lambda n:sum((-s)**j for j in range(n+1)); W=lambda n:sum((j+1)*(-s)**j for j in range(n+1))
    t=S.cancel(U(r+2)/D)
    # all-one values through 2r+1
    allb=[S.Integer(1)]
    for k in range(1,2*r+2): allb.append(S.cancel(1-s*allb[-1]-(t*allb[k-r] if k>=r else 0)))
    F=S.cancel(allb[2*r+1])
    record={'r':r,'branches':{}}
    for name,resid in [('u1',F),('u0',S.cancel(F+s))]:
        p=S.Poly(S.together(resid).as_numer_denom()[0],s,domain=S.QQ)
        # branch u0 lowers b_(2r) by 1, while b_(2r+1) is the residual.
        b=list(allb)
        if name=='u0': b[2*r]=S.cancel(b[2*r]-1); b[2*r+1]=S.cancel(F+s)
        certificates=[]
        for (a,bb),mult in p.intervals(eps=S.Rational(1,10)**16):
            if LO<a and bb<HI:
                mid=(a+bb)/2; bad=[]
                for k in range(1,2*r+1):
                    pk=S.Poly(S.together(b[k]).as_numer_denom()[0],s,domain=S.QQ)
                    if pk.count_roots(a,bb)==0 and S.sign(pk.eval(mid))<0:
                        bad.append({'k':k,'den':str(S.factor(S.together(b[k]).as_numer_denom()[1]))})
                tm=S.Poly(S.together(t-1).as_numer_denom()[0],s,domain=S.QQ)
                tn=S.Poly(S.together(t).as_numer_denom()[0],s,domain=S.QQ)
                certificates.append({'I':[str(a),str(bb)],'mult':int(mult),'negative_b':bad,
                    'tminus1_sign':str(S.sign(tm.eval(mid))) if tm.count_roots(a,bb)==0 else 'root_inside',
                    't_sign':str(S.sign(tn.eval(mid))) if tn.count_roots(a,bb)==0 else 'root_inside'})
        record['branches'][name]={
            'residual_degree':p.degree(),
            'sturm_count_closed_endpoints':int(p.count_roots(LO,HI)),
            'root_at_0':bool(p.eval(0)==0),
            'root_at_1':bool(p.eval(1)==0),
            'strict_open_root_count':len(certificates),
            'roots':certificates}
    rows.append(record)
assert Path('/Volumes/data').is_mount()
out.parent.mkdir(parents=True, exist_ok=True)
tmp=out.with_suffix('.json.tmp');tmp.write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n');os.replace(tmp,out)
print(json.dumps(rows,ensure_ascii=False,indent=2))
