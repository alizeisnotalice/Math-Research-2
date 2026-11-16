#!/usr/bin/env python3
"""Independent q=21 mask enumeration and trace replay (rules 1--7)."""
from __future__ import annotations

import argparse, copy, gzip, hashlib, json, os, sys, time
from collections import Counter
from math import gcd
from pathlib import Path

Q = 21
PS = tuple(p for p in range(2, Q) if gcd(p, Q) == 1 and Q >= 2*p)

def bad(s): raise ValueError(s)
def ilist(x,n):
    if type(x) is not list or any(type(v) is not int for v in x): bad(f'{n}: integer list required')
    return x

def k_from_mask(m): return tuple(i for i in range(2*Q) if (m >> (i % Q)) & 1)
def no_singleton(m,p):
    c=[0]*p; g=m|(m<<Q)
    for i in range(2*Q):
        if (g>>i)&1: c[i%p]+=1
    return all(x!=1 for x in c)

def source_modes(d):
    if d.get('schema')!='short-G-q-layer-input-v1' or d.get('q')!=Q or type(d.get('rows')) is not list: bad('source header')
    claimed=set(); archive=set()
    for r in d['rows']:
        if type(r) is not dict or type(r.get('p')) is not int or r.get('q')!=Q: bad('source row identity')
        p=r['p']; ss=r.get('no_mod_p_singleton_survivors')
        if p not in PS or p in claimed or type(ss) is not list: bad('source p coverage')
        claimed.add(p)
        for x in ss:
            t,k=ilist(x.get('T'),'T'),ilist(x.get('G_support'),'K')
            if t!=sorted(set(t)) or not t or t[0]!=0 or any(v<0 or v>=Q for v in t): bad('source T')
            if k!=sorted(set(t)|{v+Q for v in t}): bad('source pairing')
            key=(p,tuple(k))
            if key in archive: bad('source duplicate')
            archive.add(key)
    if claimed!=set(PS): bad('missing/extra source p')
    # Independent full bit-mask enumeration, not a DFS traversal.
    generated=set(); per_p={}
    for p in PS:
        total=eligible=survivors=0
        for m in range(1,1<<Q,2):
            total+=1; g=m|(m<<Q)
            if g & (g>>p): continue
            eligible+=1
            if no_singleton(m,p):
                generated.add((p,k_from_mask(m))); survivors+=1
        per_p[p]={'all_T_with_0':total,'p_separated_free':eligible,'no_singleton':survivors}
    if generated!=archive or len(generated)!=16070: bad('full q21 mask enumeration disagrees with input')
    return archive,per_p

def replay(row,key):
    p,Kt=key; K=set(Kt); n=max(K); t=n-Q
    if not(0<t<Q and {0,t,Q,n}<=K): bad(f'four corners {key}')
    U={0,t,Q,n}; Z=[set(),set()]; Sides=[K,{n-v for v in K}]; ec=Counter()
    for pos,e in enumerate(row.get('trace',[])):
        if type(e) is not list or not e or type(e[0]) is not str: bad(f'event shape {key}/{pos}')
        tag=e[0]
        if tag=='zero-gap':
            if len(e)!=6 or any(type(v)is not int for v in e[1:]): bad(f'gap encoding {key}/{pos}')
            _,side,i,j,h,w=e; S=Sides[side] if side in(0,1) else None
            if S is None or not(1<=i<2*Q) or i in Z[side] or j not in S or h!=i+j or not(0<=h<2*Q) or w not in S or abs(h-w)!=p: bad(f'gap witness {key}/{pos}')
            Z[side].add(i)
        elif tag=='zero-unit':
            if len(e)!=5 or any(type(v)is not int for v in e[1:]): bad(f'zero-unit encoding {key}/{pos}')
            _,side,i,j,s=e; S=Sides[side] if side in(0,1) else None; units=U if side==0 else {n-v for v in U}
            if S is None or not(1<=i<2*Q) or i in Z[side] or s not in units or j not in S or i!=s-j: bad(f'zero-unit witness {key}/{pos}')
            Z[side].add(i)
        elif tag=='unit':
            if len(e)!=3 or any(type(v)is not int for v in e[1:]): bad(f'unit encoding {key}/{pos}')
            _,side,s=e; S=Sides[side] if side in(0,1) else None; u=s if side==0 else n-s
            if S is None or s not in S or u in U or any(s-v not in Z[side] for v in S if v<s): bad(f'unit closure {key}/{pos}')
            U.add(u)
        elif tag=='unit-chain':
            if len(e)!=5 or any(type(v)is not int for v in e[1:]): bad(f'chain encoding {key}/{pos}')
            _,side,s,b,z=e; S=Sides[side] if side in(0,1) else None; units=U if side==0 else {n-v for v in U}; u=s if side==0 else n-s; ordered=sorted(S) if S else []
            if len(ordered)<3 or s!=ordered[1] or b!=ordered[2] or s in units or u in U or z not in Z[side] or not(0<z<b) or z%s: bad(f'chain witness {key}/{pos}')
            U.add(u)
        elif tag=='zero-single-product':
            if len(e)!=5 or any(type(v)is not int for v in e[1:]): bad(f'single encoding {key}/{pos}')
            _,side,s,j,k=e; S=Sides[side] if side in(0,1) else None
            # h_s=1 gives 0<=f_s<=1-g_s<1.  All other h_k terms must be Z.
            if S is None or not(0<s<2*Q) or s not in S or s in Z[side] or j not in S or k!=s+j or not(0<=k<2*Q) or any(k-v not in Z[side] for v in S if v<=k and v!=j): bad(f'single-product witness {key}/{pos}')
            Z[side].add(s)
        elif tag=='zero-two-products':
            if len(e)!=5 or any(type(v)is not int for v in e[1:]): bad(f'two encoding {key}/{pos}')
            _,side,s,t,k=e; S=Sides[side] if side in(0,1) else None
            # h_k=f_s g_t+f_t g_s, and at least one term was genuinely new.
            if S is None or not(0<s<t<2*Q) or s not in S or t not in S or k!=s+t or not(0<=k<2*Q) or (s in Z[side] and t in Z[side]) or any(k-v not in Z[side] for v in S if v<=k and v not in(s,t)): bad(f'two-product witness {key}/{pos}')
            Z[side].update((s,t))
        elif tag=='unit-single-product':
            if len(e)!=5 or any(type(v)is not int for v in e[1:]): bad(f'unit product encoding {key}/{pos}')
            _,side,s,j,k=e; S=Sides[side] if side in(0,1) else None; units=U if side==0 else {n-v for v in U}; u=s if side==0 else n-s
            if S is None or s not in S or s in units or u in U or any(s-v not in Z[side] for v in S if 0<v<s) or j not in S or k!=s+j or not(0<=k<2*Q) or any(k-v not in Z[side] for v in S if v<=k and v!=j): bad(f'unit product witness {key}/{pos}')
            U.add(u)
        else: bad(f'unknown tag {tag} at {key}/{pos}')
        ec[tag]+=1
    if ilist(row.get('unit_indices'),'U')!=sorted(U) or ilist(row.get('unknown_indices'),'unknown')!=sorted(K-U): bad(f'terminal state {key}')
    return ec

def verify(source,out):
    keys,per_p=source_modes(source)
    if out.get('schema')!='short-G-unit-propagation-v1' or type(out.get('rows')) is not list: bad('output schema')
    seen=set(); ev=Counter(); au=un=0
    for ix,r in enumerate(out['rows']):
        if type(r)is not dict or set(r)!={'p','q','G_support','unit_indices','unknown_indices','trace'}: bad(f'output row schema {ix}')
        p,q=r['p'],r['q']; K=ilist(r['G_support'],'output K')
        if type(p)is not int or q!=Q or K!=sorted(set(K)): bad(f'output identity {ix}')
        key=(p,tuple(K))
        if key not in keys or key in seen: bad(f'output coverage {ix}')
        seen.add(key); ev.update(replay(r,key)); au+=not bool(r['unknown_indices']);un+=bool(r['unknown_indices'])
    if seen!=keys or out.get('count')!=len(seen) or out.get('all_unit')!=au or out.get('unknown')!=un: bad('terminal coverage/count')
    return {'count':len(seen),'all_unit':au,'unknown':un,'event_counts':dict(ev),'mask_enumeration_per_p':per_p}

def reject(source,out,mut):
    x=copy.deepcopy(out); mut(x)
    try: verify(source,x)
    except ValueError as e:return str(e)
    raise AssertionError('semantic negative mutation accepted')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--receipt',type=Path,required=True);a=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):bad('data unmounted')
    started=time.monotonic(); sr=a.source.read_bytes();orr=a.output.read_bytes();source=json.loads(sr);out=json.loads(gzip.decompress(orr))
    if out.get('input_sha256')!=hashlib.sha256(sr).hexdigest():bad('input digest')
    result=verify(source,out)
    loc=lambda tag:next((ri,ei) for ri,r in enumerate(out['rows']) for ei,e in enumerate(r['trace']) if e[0]==tag)
    sp,tp,gp=loc('zero-single-product'),loc('zero-two-products'),loc('zero-gap')
    ur=next(i for i,r in enumerate(out['rows']) if r['unknown_indices'])
    neg={
      'single_positive_index':reject(source,out,lambda d:d['rows'][sp[0]]['trace'][sp[1]].__setitem__(2,0)),
      'single_other_term':reject(source,out,lambda d:d['rows'][sp[0]]['trace'][sp[1]].__setitem__(4,10**9)),
      'two_order':reject(source,out,lambda d:d['rows'][tp[0]]['trace'][tp[1]].__setitem__(2,d['rows'][tp[0]]['trace'][tp[1]][3])),
      'two_degree':reject(source,out,lambda d:d['rows'][tp[0]]['trace'][tp[1]].__setitem__(4,10**9)),
      'unknown_as_unit':reject(source,out,lambda d:d['rows'][ur].__setitem__('unknown_indices',[])),
      'gap_witness':reject(source,out,lambda d:d['rows'][gp[0]]['trace'][gp[1]].__setitem__(5,10**9)),
      'drop_mode':reject(source,out,lambda d:d['rows'].pop()),}
    receipt={'audit':'independent-round036-q21-mask-and-trace-v1','scope':'q=21 finite conditional necessary propagation; unknown is not feasibility and this does not settle the general conjecture.','source':str(a.source),'source_sha256':hashlib.sha256(sr).hexdigest(),'output':str(a.output),'output_sha256':hashlib.sha256(orr).hexdigest(),'audit_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'result':result,'negative_tests_rejected':neg,'interpreter':sys.version,'elapsed_seconds':time.monotonic()-started,'status':'complete'}
    a.receipt.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({'count':result['count'],'all_unit':result['all_unit'],'unknown':result['unknown'],'events':result['event_counts'],'negative':sorted(neg),'elapsed_seconds':receipt['elapsed_seconds'],'status':'complete'},sort_keys=True))
if __name__=='__main__':main()
