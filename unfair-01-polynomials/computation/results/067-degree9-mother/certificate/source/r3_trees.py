#!/usr/bin/env python3
"""Exhaust all exact r=3 annulus survivors by integer prefix trees."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

from singleton_probe import ONE, ZERO, ci, conj, extend, plus, scale


def run_tree(P,r,forward,D,depth,cap):
    rows=[((1,),(1,),(ONE,),(ONE,))]
    trace=[]
    Pb=list(map(conj,P))
    for k in range(1,depth+1):
        next_rows=[]
        for G,Q,A,B in rows:
            x=ZERO
            y=ZERO
            for j in range(1,min(9,k)+1):
                x=plus(x,scale(G[k-j],P[j]))
                y=plus(y,scale(Q[k-j],Pb[j]))
            for u,v,a,b in extend(A,B,x,y,k,D):
                if k==(9-r if forward else r) and (a!=ZERO or b!=ZERO):continue
                if forward and k==9-2*r and (u or v):continue
                next_rows.append((G+(u,),Q+(v,),A+(a,),B+(b,)))
        trace.append(len(next_rows))
        if not next_rows:return dict(status="dead",first_empty=k,layer_sizes=trace)
        if len(next_rows)>cap:return dict(status="unknown_cap",last_layer=k,layer_sizes=trace)
        rows=next_rows
    return dict(status="unknown_depth",last_layer=depth,layer_sizes=trace)


def one(job):
    row,depth,cap=job
    L,M,T,D,c,d,e,g,h=row
    P=(ONE,ci(L),ci(M),(T,1),ci(c),ci(d),ci(e),ci(g),ci(h),ONE)
    fwd=run_tree(P,3,True,D,depth,cap)
    if fwd["status"]=="dead":
        return row,dict(forward=fwd),"forward_dead"
    rev=run_tree(P[::-1],3,False,D,depth,cap)
    return row,dict(forward=fwd,reverse=rev),("reverse_dead" if rev["status"]=="dead" else "unknown")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("entry",type=Path)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=8)
    p.add_argument("--depth",type=int,default=128)
    p.add_argument("--cap",type=int,default=10000)
    a=p.parse_args()
    assert not a.out.exists()
    raw=a.entry.read_bytes()
    rows=[(row,a.depth,a.cap) for row,code in map(json.loads,raw.splitlines()) if code=="root_pass"]
    output=[]
    result=Counter()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for row,trees,status in pool.map(one,rows,chunksize=4):
            output.append((row,trees,status))
            result[status]+=1
    traces=a.out.with_suffix(".rows.jsonl")
    traces.write_text("".join(json.dumps(x,separators=(",",":"))+"\n" for x in output))
    report=dict(rows=len(rows),counts=result,depth=a.depth,cap=a.cap,
                input_sha256=hashlib.sha256(raw).hexdigest(),
                tree_sha256=hashlib.sha256(traces.read_bytes()).hexdigest())
    a.out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps(report,sort_keys=True))


if __name__=="__main__":
    main()
