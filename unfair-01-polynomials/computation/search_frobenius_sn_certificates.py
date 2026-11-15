#!/usr/bin/env python3
"""Search exact modular factorization certificates implying S_n.

For selected 0-1 C in the fixed preparation JSON, find primes with:
 irreducible; (n-1,1); and (2, odd,...,odd) with a unique 2.
All modular factors are saved. Missing patterns are status=unknown, not an
exclusion.  This script does not modify the candidate JSON or start a GPU.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as S

HERE=Path(__file__).resolve().parent
IN=HERE/'root_group_candidates_20261115.json'
OUT=HERE/'frobenius_sn_certificates_20261115.json'
x=S.symbols('x')
primes=list(S.primerange(2,1001))
names={'random_28','random_34','random_42','random_46'}

def factor_record(poly,p):
    unit, facs=S.factor_list(poly, modulus=p)
    out=[]
    for f,e in facs:
        coeff=[int(v)%p for v in reversed(S.Poly(f,x,modulus=p).all_coeffs())]
        out.append({'degree':int(S.degree(f,x)),'multiplicity':int(e),'coefficients_ascending':coeff})
    degrees=sorted([z['degree'] for z in out for _ in range(z['multiplicity'])],reverse=True)
    squarefree=all(z['multiplicity']==1 for z in out)
    return {'prime':p,'unit_mod_p':int(unit)%p,'factors':out,'cycle_type_descending':degrees,'squarefree_mod_p':squarefree}

def classify(rec,n):
    if not rec['squarefree_mod_p']: return None
    t=rec['cycle_type_descending']
    if t==[n]: return 'irreducible'
    if t==[n-1,1]: return 'nminus1_1'
    if t.count(2)==1 and all(d==2 or d%2==1 for d in t): return 'unique2_odd_rest'
    return None

def main():
    global IN, OUT
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("input",type=Path);parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args();IN,OUT=args.input,args.out
    assert Path("/Volumes/data").is_mount()
    data=json.loads(IN.read_text())
    result={'input':str(IN),'prime_limit':1000,'criterion':'irreducible + (n-1,1) + unique2/odd-rest','candidates':[]}
    for c in data['candidates']:
        if c['name'] not in names: continue
        coeff=c['coefficients_ascending']; n=c['degree']; poly=sum(a*x**i for i,a in enumerate(coeff))
        qq=S.factor_list(poly)
        rec={'name':c['name'],'degree':n,'sha256':c['sha256'],
             'coefficients_ascending':coeff,
             'qq_factor_list':[[str(f),int(e)] for f,e in qq[1]],'patterns':{},'reductions':[],'scanned_primes':0}
        if not S.Poly(poly,x).is_irreducible:
            rec['status']='unknown_reducible_over_Q';result['candidates'].append(rec)
            print(c['name'],rec['status'],[],0,flush=True)
            continue
        for p in primes:
            fr=factor_record(poly,p);rec['scanned_primes']+=1
            tag=classify(fr,n)
            if tag and tag not in rec['patterns']: rec['patterns'][tag]=fr
            if len(rec['patterns'])==3: break
        roles={'irreducible':'irreducible','nminus1_1':'long_cycle','unique2_odd_rest':'transposition_power'}
        rec['reductions']=[dict(v,role=roles[tag]) for tag,v in rec['patterns'].items()]
        rec['status']='sn_certificate_found' if len(rec['patterns'])==3 else 'unknown_missing_pattern'
        result['candidates'].append(rec)
        print(c['name'],rec['status'],sorted(rec['patterns']),rec['scanned_primes'],flush=True)
    OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__': main()
