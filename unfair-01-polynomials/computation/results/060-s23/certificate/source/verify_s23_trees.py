#!/usr/bin/env python3
"""Independently re-run the complete tree set in full integer-ring coordinates."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

from verify_s23_tail import I,Z,add,basis,integer,mulc


def canonical(row):
    return json.dumps(row,separators=(",",":")).encode()+b"\n"


def decode(rec,D):
    _,L,Ta,r,Tb,s,c,d,e,f=rec
    t=int(D%4==1)
    return [I,integer(L),((Ta-t*r)//2,r),((Tb-t*s)//2,s),
            integer(c),integer(d),integer(e),integer(f),I]


def direct(poly,D,depth,cap,remap,forced_cofactor=None):
    env=basis(D)
    choose=env["choices"]
    conj=env["conjugate"]
    cross=env["cross"]
    sign=env["sign"]
    conjugate=list(map(conj,poly))
    initial=((1,),(1,),(I,),(I,))
    states=[initial]
    trace=[]
    witnesses=[]
    for k in range(1,depth+1):
        new=[]
        failures=0
        for G,Q,A,B in states:
            def options(P,pre,slot):
                val=Z
                for i in range(1,min(k,len(P)-1)+1):
                    val=add(val,mulc(pre[k-i],P[i]))
                force=remap.get(k,(None,None))[slot]
                if force is None:
                    candidates=choose(val)
                else:
                    candidates=((force-val[0],integer(force)),) if val[1]==0 else ()
                if forced_cofactor is not None and k in forced_cofactor:
                    candidates=tuple((g,z) for g,z in candidates if g==forced_cofactor[k][slot])
                for _,z in candidates:
                    assert sign(z)>=0 and sign(add(I,mulc(-1,z)))>=0
                return candidates
            x0=cross(A,B,k)
            for g,x in options(poly,G,0):
                for q,y in options(conjugate,Q,1):
                    budget=add(add(x,y),x0)
                    if budget in (Z,I):
                        new.append((G+(g,),Q+(q,),A+(x,),B+(y,)))
                    else:
                        failures+=1
                        if len(poly)==8 and k==67:
                            witnesses.append(dict(g=g,q=q,A=x,B=y,C=budget))
        assert len({(x[0],x[1]) for x in new})==len(new)
        digest=hashlib.sha256()
        for X,Y,_,_ in sorted(new):
            digest.update(canonical((X,Y)))
        trace.append([k,len(new),failures,digest.hexdigest()])
        if not new:
            result=dict(status="dead",depth=k,trace=trace)
            if witnesses:
                result["last_witnesses"]=witnesses
            return result
        if len(new)>cap:
            return dict(status="unknown_cap",depth=k,trace=trace)
        states=new
    return dict(status="unknown_depth",depth=depth,trace=trace)


def independent(job):
    row,record,dead,depth,cap=job
    D,L,Ta,r,Tb,s,c,d,e,f=row
    P=decode(row,D)
    forward,reverse=direct(P,D,depth,cap,{5:(0,0)}),direct(P[::-1],D,depth,cap,{2:(0,0)})
    assert forward==record["forward"] and reverse==record["reverse"], row
    if dead:
        assert forward["status"]=="dead" or reverse["status"]=="dead"
        return False,None
    for [A1,B1,ea,eb],a,b in record["compatible_pairs"]:
        fa=0 if r+s*(eb-f) else None
        fb=0 if r+s*(ea-f) else None
        ga=0 if (eb-f)*(r-s*f)-s*e else None
        gb=0 if (ea-f)*(r-s*f)-s*e else None
        assert direct(P,D,depth,cap,{1:(A1,B1),5:(0,0),6:(fa,fb),7:(ga,gb)})==a
        assert direct(P[::-1],D,depth,cap,{1:(ea,eb),2:(0,0)})==b
    assert row==[3,1,4,1,4,1,0,0,0,1]
    Q=[I,Z,P[2],Z,Z,Z,Z,I]
    product=[Z]*9
    for j,x in enumerate(Q):
        product[j]=add(product[j],x)
        product[j+1]=add(product[j+1],x)
    assert product==P
    # After absorbing 1+x, the seventh-degree core has A5=B5=0 and G3=Q3=0.
    short=direct(Q,D,67,cap,{5:(0,0)},forced_cofactor={3:(0,0)})
    assert short["status"]=="dead" and short["depth"]==67
    assert short["trace"][-2][1]==1 and short["trace"][-1][2]==2
    assert {tuple(x["C"]) for x in short["last_witnesses"]}=={(2,0),(3,0)}
    return True,dict(row=row,shorter_core=Q,tree=short)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("trees",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    metadata=json.loads(a.trees.with_name("trees.json").read_bytes())
    source=a.trees.read_bytes()
    assert hashlib.sha256(source).hexdigest()==metadata["rows_sha256"]
    records=[(row,out,dead,metadata["depth_cap"],metadata["state_cap"])
             for row,out,dead in map(json.loads,source.splitlines())]
    assert len(records)==metadata["rows"]
    extras=[]
    count=Counter()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for special,witness in pool.map(independent,records,chunksize=8):
            count["dead_after_factor" if special else "dead_forward_or_reverse"]+=1
            if witness:
                extras.append(witness)
    assert len(extras)==1
    result=dict(status="verified",counts=dict(count),tree_rows=len(records),
                primary_tree_sha256=metadata["rows_sha256"],exceptional=extras[0])
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="exceptional"},sort_keys=True))


if __name__=="__main__":
    main()
