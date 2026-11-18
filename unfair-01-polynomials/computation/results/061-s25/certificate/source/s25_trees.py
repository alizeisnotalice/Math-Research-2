#!/usr/bin/env python3
"""Exact oriented prefix trees for all S25 annulus survivors."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from s25_entry import ONE,ZERO,ci,conj,plus,scale,shifts,product_cross,floor_quad


def scaled_norm(P,D):
    H=[0]*17
    for i,(u,v) in enumerate(P):
        for j,(r,s) in enumerate(P):
            H[i+j]+=u*r-v*s*D
    assert all(x%4==0 for x in H)
    return [x//4 for x in H]



def encoded(x):
    return json.dumps(x,separators=(",",":")).encode()+b"\n"


def decode(row):
    D,Ta,r,c,d,Tb,s,e,f=row
    tau=1 if D%4==1 else 2
    return [ONE,ZERO,(Ta,r*tau),ci(c),ci(d),(Tb,s*tau),ci(e),ci(f),ONE]


def known(poly,prefix,k):
    value=ZERO
    for j in range(1,min(k,8)+1):
        value=plus(value,scale(prefix[k-j],poly[j]))
    return value


def choices(base,D,forced):
    if forced is None:
        return shifts(base,D)
    if base[1] or base[0]%2:
        return ()
    return ((forced-base[0]//2,ci(forced)),)


def tree(row,P,depth,forced,cap):
    D=row[0]
    Pb=list(map(conj,P))
    H=scaled_norm(P,D)
    states={((1,),(1,)):((ONE,),(ONE,),(1,))}
    trace=[]
    for k in range(1,depth+1):
        next_layer={}
        target=forced.get(k,(None,None))
        rejected=0
        for (G,Q),(A,B,W) in states.items():
            U,V=known(P,G,k),known(Pb,Q,k)
            left,right=choices(U,D,target[0]),choices(V,D,target[1])
            cross=product_cross(A,B,k,D)
            w0=sum(G[i]*Q[k-i] for i in range(1,k))
            h0=sum(H[j]*W[k-j] for j in range(1,min(k,16)+1))
            for g,x in left:
                for q,y in right:
                    w=w0+g+q
                    assert plus(cross,plus(x,y))==ci(h0+w)
                    if h0+w not in (0,1):
                        rejected+=1
                        continue
                    key=(G+(g,),Q+(q,))
                    assert key not in next_layer
                    next_layer[key]=(A+(x,),B+(y,),W+(w,))
        digest=hashlib.sha256()
        for key in sorted(next_layer):
            digest.update(encoded(key))
        trace.append([k,len(next_layer),rejected,digest.hexdigest()])
        if not next_layer:
            return {"status":"dead","depth":k,"trace":trace}
        if len(next_layer)>cap:
            return {"status":"unknown_cap","depth":k,"trace":trace}
        states=next_layer
    return {"status":"unknown_depth","depth":depth,"trace":trace}


def one_case(job):
    row,depth,cap=job
    P=decode(row)
    results={}
    D,Ta,r,c,d,Tb,s,e,f=row
    results["forward"]=tree(row,P,depth,{1:(0,0),3:(0,0)},cap)
    g2,q2=-floor_quad(P[2],D),-floor_quad(conj(P[2]),D)
    forced={2:(0,0),4:(0 if q2 else None, 0 if g2 else None)}
    if s!=r*c:
        forced[5]=(0,0)
    results["reverse"]=tree(row,P[::-1],depth,forced,cap)
    return row,results,any(a["status"]=="dead" for a in results.values())


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("input",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    ap.add_argument("--depth",type=int,default=48)
    ap.add_argument("--cap",type=int,default=10000)
    a=ap.parse_args()
    assert not a.out.exists()
    raw=a.input.read_bytes()
    all_rows=[(tuple(row),code) for row,code in map(json.loads,raw.splitlines())]
    rows=[r for r,code in all_rows if code=="root_pass"]
    assert len(rows)==len(set(rows))
    count=Counter()
    depths=Counter()
    out=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for row,results,dead in pool.map(one_case,[(r,a.depth,a.cap) for r in rows],chunksize=8):
            count["dead" if dead else "unknown"]+=1
            for key in ("forward","reverse"):
                count[key+"_"+results[key]["status"]]+=1
                if results[key]["status"]=="dead":
                    depths[(key,results[key]["depth"])]+=1
            out.append((row,results,dead))
    payload=b"".join(map(encoded,out))
    a.out.with_suffix(".rows.jsonl").write_bytes(payload)
    result=dict(status="complete",rows=len(rows),counts=dict(count),depth_cap=a.depth,
                state_cap=a.cap,death_histogram=[[k,i,n] for (k,i),n in sorted(depths.items())],
                rows_sha256=hashlib.sha256(payload).hexdigest(),
                input_sha256=hashlib.sha256(raw).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
