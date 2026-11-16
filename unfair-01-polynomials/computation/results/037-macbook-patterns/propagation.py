#!/usr/bin/env python3
"""Necessary short-side unit propagation; no feasibility claims for unknowns."""
import argparse, gzip, hashlib, json, os, time
from pathlib import Path

def propagate(p,q,K):
    K=set(K); n=max(K); t=n-q
    U={0,t,q,n}; Z=[set(),set()]; trace=[]
    sides=[K,{n-j for j in K}]
    for side,S in enumerate(sides):
        forbidden={v for k in S for v in (k-p,k+p) if 0<=v<2*q}
        for i in range(1,2*q):
            witness=next(((j,i+j) for j in sorted(S) if i+j in forbidden),None)
            if witness:
                j,h=witness; k=next(k for k in sorted(S) if abs(k-h)==p)
                Z[side].add(i); trace.append(['zero-gap',side,i,j,h,k])
    changed=True
    while changed:
        changed=False
        for side,S in enumerate(sides):
            units=U if side==0 else {n-j for j in U}
            for s in sorted(units):
                for j in sorted(S):
                    i=s-j
                    if i>0 and i not in Z[side]:
                        Z[side].add(i);trace.append(['zero-unit',side,i,j,s]);changed=True
            for s in sorted(S-units):
                if all(s-j in Z[side] for j in S if j<s):
                    u=s if side==0 else n-s
                    U.add(u);trace.append(['unit',side,s]);changed=True
            units=U if side==0 else {n-j for j in U}
            ordered=sorted(S)
            if len(ordered)>=3:
                s,b=ordered[1:3]
                if s not in units:
                    witness=next((z for z in sorted(Z[side]) if 0<z<b and z%s==0),None)
                    if witness is not None:
                        u=s if side==0 else n-s
                        U.add(u);trace.append(['unit-chain',side,s,b,witness]);changed=True
            units=U if side==0 else {n-j for j in U}
            for s in sorted(S):
                if s==0 or s in Z[side]: continue
                for j in sorted(S):
                    k=s+j
                    if k>=2*q: continue
                    if all(k-v in Z[side] for v in S if v<=k and v!=j):
                        Z[side].add(s);trace.append(['zero-single-product',side,s,j,k]);changed=True;break
            positive=sorted(S-{0})
            for ai,s in enumerate(positive):
                for t in positive[ai+1:]:
                    k=s+t
                    if k>=2*q or (s in Z[side] and t in Z[side]): continue
                    if all(k-v in Z[side] for v in S if v<=k and v not in (s,t)):
                        Z[side].update((s,t));trace.append(['zero-two-products',side,s,t,k]);changed=True
            for s in sorted(S-units):
                # h_s=1 and every term except g_s and f_s is zero.
                if not all(s-j in Z[side] for j in S if 0<j<s):
                    continue
                for j in sorted(S):
                    k=s+j
                    if k>=2*q: continue
                    # h_k=f_s*g_j=(1-g_s)*g_j lies in [0,1).
                    if all(k-v in Z[side] for v in S if v<=k and v!=j):
                        u=s if side==0 else n-s
                        U.add(u);trace.append(['unit-single-product',side,s,j,k]);changed=True;break
    return {'p':p,'q':q,'G_support':sorted(K),'unit_indices':sorted(U),'unknown_indices':sorted(K-U),'trace':trace}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data')
    start=time.monotonic(); raw=args.input.read_bytes(); data=json.loads(gzip.decompress(raw) if args.input.suffix=='.gz' else raw)
    result=[]
    for row in data['rows']:
        for item in row['no_mod_p_singleton_survivors']:
            result.append(propagate(row['p'],row['q'],item['G_support']))
    record={'schema':'short-G-unit-propagation-v1','input_sha256':hashlib.sha256(raw).hexdigest(),'count':len(result),'all_unit':sum(not r['unknown_indices'] for r in result),'unknown':sum(bool(r['unknown_indices']) for r in result),'elapsed_seconds':time.monotonic()-start,'rows':result}
    args.output.write_bytes(gzip.compress(json.dumps(record,separators=(',',':')).encode(),mtime=0))
    print(json.dumps({k:v for k,v in record.items() if k!='rows'})); print('first unknown:',next((r for r in result if r['unknown_indices']),None)); print('p2q7:',[(r['G_support'],r['unknown_indices']) for r in result if (r['p'],r['q'])==(2,7)])
if __name__=='__main__':main()
