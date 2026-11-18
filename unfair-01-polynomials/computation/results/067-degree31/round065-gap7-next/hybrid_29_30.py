"""Hybrid projection and specialized finite-field certificates."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import sys
import time

sys.dont_write_bytecode=True
import sympy as sp
from sympy.polys import distributedmodules as dm
from sympy.polys.orderings import grevlex,lex
from prefix_29_30 import census,canon,late_codes
from modular_core_29_30 import equations,forward_eliminate
from verify_degrees_25_26 import MODS,check_modulus,pack,unpack
from verify_degrees_29_30 import exact_remaining


def do_evaluate(R,poly,z):
    result=R.zero
    c=R.domain.convert(z)
    for ex,v in poly.items():
        result+=R.from_dict({ex[:-1]+(0,):v*c**ex[-1]})
    return result


def transform(a,b,R):
    if not b:return []
    return [sum(x*column[i] for x,column in zip(a,b,strict=True)) for i in range(len(b[0]))]


def kernel_basis(R,values,gens,ordering):
    """Carry every basis element back to the polynomials supplied here."""
    if not values:
        return []
    if all(sum(sum(z) for z in f)==0 for f in values):
        k=next(i for i,v in enumerate(values) if v)
        coef=[R.one/values[k].const() if i==k else R.zero for i in range(len(values))]
        return [(R.one,coef)]
    enc=[dm.sdm_from_vector([f.as_expr()],ordering,R.domain,gens=gens) for f in values]
    B,M=dm.sdm_groebner(enc,dm.sdm_nf_mora,ordering,R.domain,extended=True)
    result=[]
    for b,m in zip(B,M,strict=True):
        f=R.from_expr(dm.sdm_to_vector(b,gens,R.domain,n=1)[0])
        h=[R.from_expr(q) for q in dm.sdm_to_vector(m,gens,R.domain,n=len(values))]
        assert f==sum(x*y for x,y in zip(h,values,strict=True))
        result.append((f,h))
    return result


def prelude(j,shape,mod):
    R,vs,E=equations(j,shape,mod)
    Q=forward_eliminate(R,vs,E,True)
    gens=tuple(v.as_expr() for v in vs)
    poly=[f for f,h in Q]
    multipliers=[h for f,h in Q]
    used={i for i in range(len(vs)-1) if any(any(g[i] for g in f) for f in poly)}
    projections=[]
    if not used:
        projections=Q
    elif len(used)<=2:
        for f,m in kernel_basis(R,poly,gens,lex):
            if all(not any(q[:-1]) for q in f):
                projections.append((f,transform(m,multipliers,R)))
    for f,v in projections:
        assert f==sum(x*y for x,y in zip(v,E,strict=True))
    return R,vs,E,Q,projections


def specific(data,z):
    R,vs,E,Q,_=data
    arr=[]
    for f,h in Q:
        j=do_evaluate(R,f,z)
        if j:
            arr.append((j,[do_evaluate(R,x,z) for x in h]))
    base=[f for f,h in arr]
    pivots=forward_eliminate(R,vs,base,True)
    L=[g for g,h in pivots]
    # Directly recognize a constant before invoking the solver.
    for i,f in enumerate(L):
        if f.is_ground:
            relative=[q/f.const() for q in pivots[i][1]]
            orig=transform(relative,[h for g,h in arr],R)
            return orig
    used={i for i in range(len(vs)-1) if any(any(g[i] for g in f) for f in L)}
    if not used:return None
    ordered=tuple(vs[i].as_expr() for i in sorted(used))
    B=kernel_basis(R,L,ordered,grevlex)
    for f,h in B:
        if f.is_ground:
            rel=transform([u/f.const() for u in h],[y for x,y in pivots],R)
            return transform(rel,[h for f,h in arr],R)
    return None


def verify_projected(R,E,proof_list):
    data=[]
    for row in proof_list:
        f=unpack(row['f'],R)
        assert f and all(not any(x[:-1]) for x in f)
        vs=[unpack(h,R) for h in row['h']]
        assert len(E)==len(vs)
        assert f==sum(x*y for x,y in zip(vs,E,strict=True))
        data.append(f)
    return data


def verify_unit(R,E,z,proof):
    h=[unpack(f,R) for f in proof]
    assert len(h)==len(E)
    assert sum(x*do_evaluate(R,y,z) for x,y in zip(h,E,strict=True))==R.one


def modvalue(word,mod):
    p,_,weight=mod
    return sum(x*y for x,y in zip(weight,word,strict=True))%p


def filter_field(j,s,mod,words,existing=None):
    p=mod[0]
    if existing is None:
        stage=prelude(j,s,mod)
        R,vs,E,_,projs=stage
        kept={'projections':[{'f':pack(f,p),'h':[pack(x,p) for x in h]} for f,h in projs],
              'units':{}}
    else:
        R,vs,E=equations(j,s,mod)
        kept=existing
    f=verify_projected(R,E,kept['projections'])
    live=[]
    tried={}
    projected=specific_count=0
    seen_units=set()
    for w in words:
        z=modvalue(w,mod)
        if any(do_evaluate(R,g,z) for g in f):
            projected+=1
            continue
        key=str(z)
        if key not in tried:
            if existing is None:
                v=specific(stage,z)
                if v is not None:
                    kept['units'][key]=[pack(x,p) for x in v]
            tried[key]=kept['units'].get(key)
            if tried[key] is not None:
                verify_unit(R,E,z,tried[key])
        if tried[key] is None:
            live.append(w)
        else:
            specific_count+=1
            seen_units.add(key)
    assert seen_units==set(kept['units'])
    counts=dict(projected_words=projected,specialized_words=specific_count,projected_proofs=len(f),
                specialized_proofs=len(seen_units),live_words=len(live))
    return live,kept,counts


def run_degree(j,base,already=None,offset=0,limit=None):
    saved=None if already is None else {canon(x['shape']):x['fields'] for x in already}
    counts=Counter()
    retained=[]
    certificates=[]
    rows=base['late'][offset:] if limit is None else base['late'][offset:offset+limit]
    for idx,r in enumerate(rows,offset+1):
        shape=canon(r['shape'])
        words=late_codes(j,shape,r['A'],r['B'],counts)
        counts['shapes']+=1
        counts['words']+=len(words)
        result={}
        for mod in MODS:
            if not words:
                break
            key=str(mod[0])
            before=len(words)
            words,witness,after=filter_field(j,shape,mod,words,None if saved is None else saved[shape][key])
            result[key]=witness
            counts['mod'+key+'_rejected']+=before-len(words)
            for tag,count in after.items():
                counts['mod'+key+'_'+tag]+=count
        if words:
            print(j,idx,'unresolved',shape,words,flush=True)
            retained.append(exact_remaining(j,shape,words))
            counts['unresolved']+=len(words)
        certificates.append({'shape':shape,'fields':result,'surviving_words':words})
        if idx%300==0:
            print(j,idx,'of',len(base['late']),dict(counts),flush=True)
    assert counts['words']+counts['tail_interval_rejected']==counts['support_word_count']
    assert sum(counts[f'mod{p}_rejected'] for p,*_ in MODS)+counts['unresolved']==counts['words']
    return {'degree':j,'early_counts':base['counts'],'late_counts':dict(counts),'certificates':certificates,'unresolved':retained}


def main():
    if not __debug__:raise SystemExit('requires assertions')
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--only',type=int,choices=(29,30),required=True)
    ap.add_argument('--cache',type=Path)
    ap.add_argument('--witness',type=Path)
    ap.add_argument('--offset',type=int,default=0)
    ap.add_argument('--limit',type=int)
    args=ap.parse_args()
    assert not (args.cache and args.witness) and not args.out.exists()
    assert args.offset>=0 and (args.limit is None or args.limit>0)
    start=time.monotonic()
    j=args.only
    for p,t,l in MODS:check_modulus(p,t,l)
    b=census(j) if args.cache is None else next(x for x in json.loads(args.cache.read_text()) if x['degree']==j)
    saved=None if not args.witness else json.loads(args.witness.read_text())['result']['certificates']
    data=run_degree(j,b,saved,args.offset,args.limit)
    here=Path(__file__).parent
    files=('verify_degrees_17_18.py','verify_degree_19.py','verify_degrees_21_22.py','verify_degrees_23_24.py',
           'verify_degrees_25_26.py','prefix_29_30.py','modular_core_29_30.py','verify_degrees_29_30.py',Path(__file__).name)
    r={'status':'complete','mode':'verify' if args.witness else 'generate','result':data,'offset':args.offset,
       'limit':args.limit,'read_early_cache':args.cache is not None,'seconds':time.monotonic()-start,
       'python':sys.version,'sympy':sp.__version__,
       'sources_sha256':{n:sha256((here/n).read_bytes()).hexdigest() for n in files}}
    if args.witness:r['witness_sha256']=sha256(args.witness.read_bytes()).hexdigest()
    with args.out.open('x') as f:f.write(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
    print(j,data['late_counts'],r['seconds'],flush=True)


if __name__=='__main__':main()
