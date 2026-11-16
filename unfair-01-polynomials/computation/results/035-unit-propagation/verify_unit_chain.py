#!/usr/bin/env python3
"""Independent replay of the enhanced round-35 unit-chain certificate.

Standalone by design: it does not import the producer or the earlier audit.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import os
import sys
import time
from collections import Counter
from math import gcd
from pathlib import Path


def bad(s: str):
    raise ValueError(s)


def ilist(x, name):
    if type(x) is not list or any(type(v) is not int for v in x):
        bad(f"{name} must be a list of exact integers")
    return x


def source_keys(data):
    if data.get("schema") != "short-G-singleton-support-v1": bad("source schema")
    if data.get("range") != {"q<=":20,"coprime_1<p<q":True,"q_ge_2p_only":False}: bad("source range")
    pairs=set(); keys=set()
    for r in data.get("rows",[]):
        p,q=r.get("p"),r.get("q")
        if type(p) is not int or type(q) is not int or not (1<p<q and gcd(p,q)==1) or (p,q) in pairs: bad("source p/q")
        pairs.add((p,q)); survivors=r.get("no_mod_p_singleton_survivors")
        if type(survivors) is not list or r.get("no_mod_p_singleton_count") != len(survivors): bad("source count")
        for x in survivors:
            t,k=ilist(x.get("T"),"source T"),ilist(x.get("G_support"),"source K")
            if t != sorted(set(t)) or not t or t[0]!=0 or any(v<0 or v>=q for v in t): bad("source T shape")
            if k != sorted(set(t)|{v+q for v in t}): bad("source pairing")
            if any(v+p in set(k) for v in k): bad("source p gap")
            if any(c==1 for c in Counter(v%p for v in k).values()): bad("source singleton")
            key=(p,q,tuple(k))
            if key in keys: bad("source duplicate")
            keys.add(key)
    target={(p,q) for q in range(3,21) for p in range(2,q) if gcd(p,q)==1}
    if pairs != target or len(keys)!=18734 or data.get("summary",{}).get("no_singleton_survivors")!=18734: bad("source coverage")
    return keys


def replay(row,key):
    p,q,kt=key; K=set(kt); n=max(K); t=n-q
    if not (0<t<q and {0,t,q,n}<=K): bad(f"four corners {key}")
    U={0,t,q,n}; Z=[set(),set()]; Sides=[K,{n-v for v in K}]; ec=Counter()
    for pos,e in enumerate(row.get("trace",[])):
        if type(e) is not list or not e or type(e[0]) is not str: bad(f"event shape {key} {pos}")
        tag=e[0]
        if tag=='zero-gap':
            if len(e)!=6 or any(type(x) is not int for x in e[1:]): bad(f"gap encoding {key} {pos}")
            _,side,i,j,h,w=e; S=Sides[side] if side in (0,1) else None
            if S is None or not(1<=i<2*q) or i in Z[side] or j not in S or h!=i+j or not(0<=h<2*q) or w not in S or abs(h-w)!=p: bad(f"gap witness {key} {pos}")
            Z[side].add(i)
        elif tag=='zero-unit':
            if len(e)!=5 or any(type(x) is not int for x in e[1:]): bad(f"zero-unit encoding {key} {pos}")
            _,side,i,j,s=e; S=Sides[side] if side in (0,1) else None; units=U if side==0 else {n-v for v in U}
            if S is None or not(1<=i<2*q) or i in Z[side] or s not in units or j not in S or i!=s-j: bad(f"zero-unit witness {key} {pos}")
            Z[side].add(i)
        elif tag=='unit':
            if len(e)!=3 or any(type(x) is not int for x in e[1:]): bad(f"unit encoding {key} {pos}")
            _,side,s=e; S=Sides[side] if side in (0,1) else None; u=s if side==0 else n-s
            if S is None or s not in S or u in U or any(s-v not in Z[side] for v in S if v<s): bad(f"unit closure {key} {pos}")
            U.add(u)
        elif tag=='unit-single-product':
            if len(e)!=5 or any(type(x) is not int for x in e[1:]): bad(f"product encoding {key} {pos}")
            _,side,s,j,k=e; S=Sides[side] if side in (0,1) else None; units=U if side==0 else {n-v for v in U}; u=s if side==0 else n-s
            if (S is None or s not in S or s in units or u in U or any(s-v not in Z[side] for v in S if 0<v<s)
                or j not in S or k!=s+j or not(0<=k<2*q) or any(k-v not in Z[side] for v in S if v<=k and v!=j)):
                bad(f"single-product witness {key} {pos}")
            U.add(u)
        elif tag=='unit-chain':
            if len(e)!=5 or any(type(x) is not int for x in e[1:]): bad(f"chain encoding {key} {pos}")
            _,side,s,b,z=e; S=Sides[side] if side in (0,1) else None; units=U if side==0 else {n-v for v in U}; u=s if side==0 else n-s
            ordered=sorted(S) if S is not None else []
            # s,b are the first two strictly positive G-support indices. Thus
            # below b only 0,s occur in G, exactly what the recurrence needs.
            if (len(ordered)<3 or s!=ordered[1] or b!=ordered[2] or s in units or u in U
                or z not in Z[side] or not(0<z<b) or z % s != 0):
                bad(f"chain witness {key} {pos}")
            U.add(u)
        else: bad(f"unknown tag {tag} {key} {pos}")
        ec[tag]+=1
    if ilist(row.get("unit_indices"),"U") != sorted(U) or ilist(row.get("unknown_indices"),"unknown") != sorted(K-U): bad(f"terminal state {key}")
    return ec


def verify(source, out):
    keys=source_keys(source)
    if out.get("schema")!='short-G-unit-propagation-v1' or type(out.get("rows")) is not list: bad("output schema")
    seen=set(); events=Counter(); au=unk=0
    for pos,row in enumerate(out['rows']):
        if type(row) is not dict or set(row)!={'p','q','G_support','unit_indices','unknown_indices','trace'}: bad(f"row schema {pos}")
        p,q=row['p'],row['q']; k=ilist(row['G_support'],'output K')
        if type(p) is not int or type(q) is not int or k!=sorted(set(k)): bad(f"identity {pos}")
        key=(p,q,tuple(k))
        if key not in keys or key in seen: bad(f"coverage item {pos}")
        seen.add(key); events.update(replay(row,key)); au += not bool(row['unknown_indices']); unk += bool(row['unknown_indices'])
    if seen!=keys or out.get('count')!=len(seen) or out.get('all_unit')!=au or out.get('unknown')!=unk: bad("final coverage/count")
    return {'count':len(seen),'all_unit':au,'unknown':unk,'event_counts':dict(events)}


def reject(source,out,mutate):
    x=copy.deepcopy(out); mutate(x)
    try: verify(source,x)
    except ValueError as e: return str(e)
    raise AssertionError('negative mutation accepted')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--source',type=Path,required=True);ap.add_argument('--chain',type=Path,required=True);ap.add_argument('--receipt',type=Path,required=True);a=ap.parse_args()
    if not os.path.ismount('/Volumes/data'): bad('data unmounted')
    started=time.monotonic(); sr=a.source.read_bytes(); cr=a.chain.read_bytes(); source=json.loads(sr); out=json.loads(gzip.decompress(cr))
    if out.get('input_sha256')!=hashlib.sha256(sr).hexdigest(): bad('input digest binding')
    result=verify(source,out)
    chainloc=next((ri,ei) for ri,r in enumerate(out['rows']) for ei,e in enumerate(r['trace']) if e[0]=='unit-chain')
    gaploc=next((ri,ei) for ri,r in enumerate(out['rows']) for ei,e in enumerate(r['trace']) if e[0]=='zero-gap')
    neg={
      'chain_nonminimal_s':reject(source,out,lambda d:d['rows'][chainloc[0]]['trace'][chainloc[1]].__setitem__(2,999)),
      # The first saved chain event has (s,b,z)=(3,7,6), and its prior
      # zero-gap trace has also established z=5.  Thus 5 remains an actual
      # zero in the required interval but is genuinely not divisible by s.
      'chain_zero_not_divisible':reject(source,out,lambda d:d['rows'][chainloc[0]]['trace'][chainloc[1]].__setitem__(4,5)),
      # z=3 stays in 0<z<b and is divisible by s, but it has not entered
      # side-0 Z at this point (it only appears as a later zero-unit result).
      'chain_missing_zero':reject(source,out,lambda d:d['rows'][chainloc[0]]['trace'][chainloc[1]].__setitem__(4,3)),
      'gap_witness':reject(source,out,lambda d:d['rows'][gaploc[0]]['trace'][gaploc[1]].__setitem__(5,999)),
      'drop_row':reject(source,out,lambda d:d['rows'].pop()),
    }
    receipt={'audit':'independent-round035-unit-chain-replay-v1','scope':'Conditional necessary propagation only; no factorization or conjecture conclusion.',
      'source':str(a.source),'source_sha256':hashlib.sha256(sr).hexdigest(),'chain':str(a.chain),'chain_sha256':hashlib.sha256(cr).hexdigest(),
      'audit_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'result':result,'negative_tests_rejected':neg,
      'interpreter':sys.version,'elapsed_seconds':time.monotonic()-started,'status':'complete'}
    a.receipt.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n'); print(json.dumps({'result':result,'negative_tests':sorted(neg),'elapsed_seconds':receipt['elapsed_seconds'],'status':'complete'},sort_keys=True))


if __name__=='__main__': main()
