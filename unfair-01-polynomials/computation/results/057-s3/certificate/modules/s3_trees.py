"""Exact prefix trees for degree-eight mothers with nonrational position 3."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import hashlib
import json
import math
import time

ZERO = (0, 0)
ONE = (1, 0)


def plus(x, y):
    return x[0]+y[0], x[1]+y[1]


def times(x, y, trace, norm):
    a,b=x
    c,d=y
    return a*c-b*d*norm, a*d+b*c+b*d*trace


def sign(x, trace, discriminant):
    a,b=x
    u=2*a+b*trace
    if b==0:
        return (u>0)-(u<0)
    if b>0 and u>=0:
        return 1
    if b<0 and u<=0:
        return -1
    delta=u*u-b*b*discriminant
    return ((delta>0)-(delta<0))*(1 if u>0 else -1)


def shifts(base, trace, discriminant, row, reverse):
    a,b=base
    if row == (3 if reverse else 5):
        return [(-a,ZERO)] if b==0 else []
    if row==2 and not reverse:
        if sign(base,trace,discriminant)>=0 and sign((1-a,-b),trace,discriminant)>=0:
            return [(0,base)]
        return []
    if b==0:
        return [(-a,ZERO),(1-a,ONE)]
    s=math.isqrt(b*b*discriminant)
    if b<0:
        s=-s-1
    h=-((2*a+b*trace+s)//2)
    output=(a+h,b)
    assert sign(output,trace,discriminant)>0 and sign((1-output[0],-b),trace,discriminant)>0
    return [(h,output)]


def setup(row, reverse):
    L,M,T,N,D,b,c,d,e=row
    assert M==L*L and T*T-4*N==D>0 and math.isqrt(D)**2!=D
    P=[ONE,(L,0),(M,0),(0,1),(b,0),(c,0),(d,0),(e,0),ONE]
    Pb=[(a+T*z,-z) for a,z in P]
    if reverse:
        P.reverse();Pb.reverse()
    return P,Pb,T,N,D


def cofactor_base(P, G, k):
    terms=[(G[k-i],P[i]) for i in range(1,min(k,8)+1)]
    return sum(q*p for q,(p,_) in terms),sum(q*z for q,(_,z) in terms)


def prefix_hash(states):
    encoded=json.dumps(sorted((s[0],s[1]) for s in states),separators=(',',':')).encode()
    return hashlib.sha256(encoded).hexdigest()


def traverse(row, reverse, algorithm, max_depth, cap):
    P,Pb,T,N,D=setup(row,reverse)
    if algorithm=='ab':
        states=[((1,),(1,),(ONE,),(ONE,))]
        H=None
    else:
        states=[((1,),(1,),(1,))]
        # Build H as F^2+T*x^r*F+N*x^(2r), independently of the AB multiplication.
        position=5 if reverse else 3
        f=[p for p,_ in P]
        h=[0]*17
        for i,x in enumerate(f):
            for j,y in enumerate(f):h[i+j]+=x*y
            h[i+position]+=T*x
        h[2*position]+=N
        H=h
    summaries=[]
    trace=[]
    for k in range(1,max_depth+1):
        out=[]
        rejected=0
        for state in states:
            G,Q=state[:2]
            u=cofactor_base(P,G,k)
            v=cofactor_base(Pb,Q,k)
            if algorithm=='ab':
                A,B=state[2:]
                cross=ZERO
                for i in range(1,k):
                    cross=plus(cross,times(A[i],B[k-i],T,N))
            else:
                W=state[2]
                subtotal=sum(G[i]*Q[k-i] for i in range(1,k))
                cross= subtotal+sum(H[i]*W[k-i] for i in range(1,min(16,k)+1))
            for g,ak in shifts(u,T,D,k,reverse):
                for q,bk in shifts(v,T,D,k,reverse):
                    if algorithm=='ab':
                        if plus(cross,plus(ak,bk)) not in (ZERO,ONE):
                            rejected+=1;continue
                        out.append((G+(g,),Q+(q,),A+(ak,),B+(bk,)))
                    else:
                        if cross+g+q not in (0,1):
                            rejected+=1;continue
                        out.append((G+(g,),Q+(q,),W+(subtotal+g+q,)))
        assert len({(s[0],s[1]) for s in out})==len(out)
        trace.append((k,len(out),rejected,prefix_hash(out)))
        summaries.append(len(out))
        if not out:
            return 'empty',trace,summaries
        if len(out)>cap:
            return 'unresolved_state_cap',trace,summaries
        states=out
    return 'unresolved_depth',trace,summaries


def check(job):
    row,depth,cap=job
    layers=0
    transcripts=[]
    for reverse in (False,True):
        primary=traverse(row,reverse,'ab',depth,cap)
        independent=traverse(row,reverse,'hgq',depth,cap)
        assert primary==independent,(row,reverse)
        status,trace,counts=primary
        transcripts.append(trace)
        layers+=len(trace)
        if status=='empty':
            return row,'empty_reverse' if reverse else 'empty_forward',layers,len(trace),max(counts),hashlib.sha256(json.dumps(transcripts,separators=(',',':')).encode()).hexdigest()
    return row,'unresolved',layers,0,0,hashlib.sha256(json.dumps(transcripts,separators=(',',':')).encode()).hexdigest()


def main():
    if not __debug__:
        raise RuntimeError('run without -O')
    ap=ArgumentParser()
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--depth',type=int,default=96)
    ap.add_argument('--state-cap',type=int,default=10000)
    a=ap.parse_args()
    assert not a.out.exists() and 1<=a.workers<=13 and a.depth>=1 and a.state_cap>=1
    a.out.mkdir()
    raw=a.input.read_bytes()
    rows=[tuple(json.loads(line)) for line in raw.splitlines()]
    assert rows and len(rows)==len(set(rows)) and all(len(row)==9 and all(type(v) is int for v in row) for row in rows)
    countdown=Counter()
    deaths=Counter()
    retained=[]
    start=time.monotonic()
    total_layers=0
    with (a.out/'rows.jsonl').open('xb') as ledger, ProcessPoolExecutor(a.workers) as pool:
        for row,status,layers,death,live,digest in pool.map(check,((row,a.depth,a.state_cap) for row in rows),chunksize=128):
            countdown[status]+=1
            total_layers+=layers
            if status=='unresolved':retained.append(row)
            else:deaths[death]+=1
            ledger.write(json.dumps([row,status,layers,death,live,digest],separators=(',',':')).encode()+b'\n')
    assert sum(countdown.values())==len(rows) and sum(deaths.values())+len(retained)==len(rows)
    result={'status':'complete_necessary_trees','input_sha256':hashlib.sha256(raw).hexdigest(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'ledger_sha256':hashlib.sha256((a.out/'rows.jsonl').read_bytes()).hexdigest(),
            'mothers':len(rows),'states':dict(countdown),'deaths':dict(deaths),'layers':total_layers,
            'unresolved':retained,'seconds':time.monotonic()-start,'workers':a.workers,'depth':a.depth,'state_cap':a.state_cap}
    (a.out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('mothers','states','deaths','layers','seconds')}))


if __name__=='__main__':
    main()
