#!/usr/bin/env python3
"""Independently re-run the complete tree set in full integer-ring coordinates."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

from verify_s34_entry import I,Z,add,basis,integer,mulc


def canonical(row):
    return json.dumps(row,separators=(",",":")).encode()+b"\n"


def decode(rec,D):
    _,L,M,Ta,r,u,c,d,e=rec
    t=int(D%4==1)
    a=((Ta-t*r)//2,r)
    b=add(mulc(L,a),integer(u))
    return [I,integer(L),integer(M),a,b,integer(c),integer(d),integer(e),I]


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
    D,L,M,Ta,r,u,c,d,e=row
    P=decode(row,D)
    fwd={1:(0,0),4:(0,0)}
    if abs(L)>=2:
        fwd.update({2:(0,0),5:(0,0)})
    env=basis(D)
    choices,conj=env["choices"],env["conjugate"]
    left=choices(P[3])[0][0]+M*L
    right=choices(conj(P[3]))[0][0]+M*L
    forced={3:(0,0),6:(0 if right else None,0 if left else None)}
    if M not in (0,1):
        forced.update({1:(0,0),5:(0,0)})
    forward=direct(P,D,depth,cap,fwd)
    reverse=direct(P[::-1],D,depth,cap,forced)
    assert dead and forward["status"]=="dead"
    assert forward==record["forward"] and reverse==record["reverse"], row
    rev=reverse["depth"] if reverse["status"]=="dead" else -1
    return forward["depth"],rev


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
    fwd,rev=Counter(),Counter()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for x,y in pool.map(independent,records,chunksize=2):
            fwd[x]+=1
            rev[y]+=1
    assert sum(fwd.values())==sum(rev.values())==len(records)
    result=dict(status="verified",rows=len(records),primary_tree_sha256=metadata["rows_sha256"],
                forward=dict(sorted(fwd.items())),reverse=dict(sorted(rev.items())))
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
