#!/usr/bin/env python3
"""Exact necessary rows 2--4 for degree-eight S={1} cores.

No parent tree/search code is imported.  Elements are p+q*a with
a^2=T*a-N and a=(T+sqrt(Delta))/2; all sign/floor decisions use integers.
The result is a necessary-prefix ledger only: no survivor is a candidate.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

Pair = tuple[int, int]  # p+q*a
ZERO: Pair = (0, 0)
ONE: Pair = (1, 0)

BOUNDS = {2: (17, 1), 3: (65, 2), 4: (257, 4)}


def add(x: Pair, y: Pair) -> Pair:
    return (x[0] + y[0], x[1] + y[1])


def neg(x: Pair) -> Pair:
    return (-x[0], -x[1])


def scale(k: int, x: Pair) -> Pair:
    return (k*x[0], k*x[1])


def mul(x: Pair, y: Pair, T: int, N: int) -> Pair:
    p,q = x; r,s = y
    return (p*r-q*s*N, p*s+q*r+q*s*T)


def conj(x: Pair, T: int) -> Pair:
    p,q = x
    return (p+q*T, -q)


def sign(x: Pair, T: int, Delta: int) -> int:
    """Exact sign of p+q*(T+sqrt(Delta))/2."""
    p,q = x
    r = 2*p + q*T
    if q == 0:
        return (r > 0) - (r < 0)
    if r >= 0 and q > 0:
        return 1
    if r <= 0 and q < 0:
        return -1
    left = r*r
    right = q*q*Delta
    if r > 0:  # q<0
        return (left > right) - (left < right)
    # r<0, q>0
    return (right > left) - (right < left)


def floor_pair(x: Pair, T: int, Delta: int) -> int:
    p,q = x
    B = abs(p) + abs(q)*(abs(T)+math.isqrt(Delta)+2) + 3
    lo,hi = -B,B
    while hi-lo > 1:
        m = (lo+hi)//2
        if sign((p-m,q),T,Delta) >= 0:
            lo = m
        else:
            hi = m
    return lo


def unit_choices(base: Pair, T: int, Delta: int) -> list[tuple[int, Pair]]:
    """All integer h with 0 <= base+h <= 1, retaining rational endpoints."""
    p,q = base
    if q == 0:
        return [(-p,(0,0)), (1-p,(1,0))]
    h = -floor_pair(base,T,Delta)
    z = (p+h,q)
    assert sign(z,T,Delta) > 0 and sign(add((1,0),neg(z)),T,Delta) > 0
    return [(h,z)]


def within(x: Pair, num: int, den: int, T: int, Delta: int) -> bool:
    """Strict |x| < num/den under the chosen real embedding."""
    p,q = x
    return (sign((num-den*p,-den*q),T,Delta) > 0 and
            sign((num+den*p, den*q),T,Delta) > 0)


def both_within(x: Pair, num: int, den: int, T: int, Delta: int) -> bool:
    return within(x,num,den,T,Delta) and within(conj(x,T),num,den,T,Delta)


def c_binary(A: tuple[Pair,...], B: tuple[Pair,...], k: int, T: int, N: int) -> int | None:
    z = ZERO
    for i in range(k+1):
        z = add(z,mul(A[i],B[k-i],T,N))
    return z[0] if z[1] == 0 and z[0] in (0,1) else None


def powers_for_newton(a: Pair, b: int, c: int, d: int, T: int, N: int) -> tuple[Pair,Pair,Pair,Pair]:
    """Newton sums through s4 for P=1+a*x+b*x2+c*x3+d*x4+... ."""
    s1 = neg(a)
    s2 = add(neg(mul(a,s1,T,N)),(-2*b,0))
    s3 = add(add(neg(mul(a,s2,T,N)),neg(scale(b,s1))),(-3*c,0))
    s4 = add(add(add(neg(mul(a,s3,T,N)),neg(scale(b,s2))),neg(scale(c,s1))),(-4*d,0))
    return s1,s2,s3,s4


def core_list() -> list[tuple[int,int,int]]:
    out=[]
    for T in range(-19,20):
        for N in range(-99,100):
            D=T*T-4*N
            if D<=0 or math.isqrt(D)**2 == D: continue
            if 20-abs(T)<=0 or D >= (20-abs(T))**2: continue
            out.append((T,N,D))
    assert len(out) == 1140
    return out


def atomic_json(path: Path, obj: object) -> None:
    tmp=path.with_suffix(path.suffix+'.partial')
    with tmp.open('w') as f:
        json.dump(obj,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one(core: tuple[int,int,int], outdir: str, code_hash: str, max_seconds: float) -> dict:
    T,N,Delta=core; root=Path(outdir); rec=root/'cores'/f'T{T}_N{N}.json'
    binding={'code_sha256':code_hash,'core':[T,N,Delta],'bounds':BOUNDS}
    if rec.exists():
        old=json.loads(rec.read_text())
        assert old['binding']==binding and old['status']=='complete'
        return old['summary'] | {'cached':True}
    began=time.monotonic(); a=(0,1); abar=conj(a,T)
    g1=-floor_pair(a,T,Delta); q1=-floor_pair(abar,T,Delta)
    A1=(g1,1); B1=(q1+T,-1)
    assert sign(A1,T,Delta)>0 and sign(add(ONE,neg(A1)),T,Delta)>0
    assert sign(B1,T,Delta)>0 and sign(add(ONE,neg(B1)),T,Delta)>0
    assert c_binary((ONE,A1),(ONE,B1),1,T,N)==1
    counts=Counter(cores=1,b_box=83); states=[]
    # Newton s2 and complete row 2.
    for b in range(-41,42):
        _,s2,_,_=powers_for_newton(a,b,0,0,T,N)
        if not both_within(s2,*BOUNDS[2],T,Delta): counts['s2_reject']+=1;continue
        counts['s2_pass']+=1
        for rho,A2 in unit_choices(scale(g1,a),T,Delta):
            for sigma,B2 in unit_choices(scale(q1,abar),T,Delta):
                if c_binary((ONE,A1,A2),(ONE,B1,B2),2,T,N) is None: counts['c2_reject']+=1;continue
                counts['row2_states']+=1;states.append((b,rho,sigma,A2,B2))
    # Newton s3, row 3, then Newton s4 and row 4.
    final=[]
    for b,rho,sigma,A2,B2 in states:
        for c in range(-92,93):
            _,_,s3,_=powers_for_newton(a,b,c,0,T,N)
            if not both_within(s3,*BOUNDS[3],T,Delta): counts['s3_reject']+=1;continue
            counts['s3_pass']+=1
            U3=add(scale(rho,a),scale(b,add((g1,0),neg(a))))
            V3=add(scale(sigma,abar),scale(b,add((q1,0),neg(abar))))
            for eta,A3 in unit_choices(U3,T,Delta):
                for zeta,B3 in unit_choices(V3,T,Delta):
                    if c_binary((ONE,A1,A2,A3),(ONE,B1,B2,B3),3,T,N) is None: counts['c3_reject']+=1;continue
                    counts['row3_states']+=1
                    for d in range(-120,121):
                        _,_,_,s4=powers_for_newton(a,b,c,d,T,N)
                        if not both_within(s4,*BOUNDS[4],T,Delta): counts['s4_reject']+=1;continue
                        counts['s4_pass']+=1
                        U4=add(scale(eta,a),scale(c,add((g1,0),neg(a))))
                        V4=add(scale(zeta,abar),scale(c,add((q1,0),neg(abar))))
                        for alpha,A4 in unit_choices(U4,T,Delta):
                            for beta,B4 in unit_choices(V4,T,Delta):
                                if c_binary((ONE,A1,A2,A3,A4),(ONE,B1,B2,B3,B4),4,T,N) is None: counts['c4_reject']+=1;continue
                                counts['row4_states']+=1
                                final.append([b,c,d,rho,sigma,eta,zeta,alpha,beta])
        if time.monotonic()-began > max_seconds:
            raise TimeoutError(f'core time cap T={T},N={N}')
    assert len(final)==counts['row4_states']
    summary={'T':T,'N':N,'Delta':Delta,'counts':dict(counts),'row4_states':len(final),'seconds':time.monotonic()-began}
    atomic_json(rec,{'status':'complete','binding':binding,'summary':summary,'records':final})
    return summary


def symbolic_check() -> None:
    import sympy as sp
    a,b,c,d=sp.symbols('a b c d'); g1,rho,eta,alpha=sp.symbols('g1 rho eta alpha')
    G2=rho-b; G3=eta-c; G4=alpha-d-b*rho+b*b
    A2=sp.expand(G2+a*g1+b)
    A3=sp.expand(G3+a*G2+b*g1+c)
    A4=sp.expand(G4+a*G3+b*G2+c*g1+d)
    assert sp.expand(A2-(a*g1+rho)) == 0
    assert sp.expand(A3-(a*rho+b*(g1-a)+eta)) == 0
    assert sp.expand(A4-(a*eta+c*(g1-a)+alpha)) == 0


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--core-limit',type=int,default=None)
    ap.add_argument('--core-seconds',type=float,default=60.0)
    ap.add_argument('--seconds',type=float,default=14400.0)
    ap.add_argument('--symbolic-check',action='store_true')
    args=ap.parse_args()
    if args.symbolic_check:
        symbolic_check();print(json.dumps({'status':'symbolic_formulas_passed'}));return
    if not __debug__: raise RuntimeError('refuse python -O')
    assert 1<=args.workers<=4
    if args.out.exists() and not (args.out/'cores').is_dir(): raise FileExistsError(args.out)
    args.out.mkdir(parents=True,exist_ok=True);(args.out/'cores').mkdir(exist_ok=True)
    code_hash=sha(Path(__file__));cores=core_list();
    core_list_sha256=hashlib.sha256(json.dumps(cores,separators=(',',':')).encode()).hexdigest()
    if args.core_limit is not None: cores=cores[:args.core_limit]
    start=time.monotonic();done=[];total=Counter();deadline=start+args.seconds
    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs=[pool.submit(one,x,str(args.out),code_hash,args.core_seconds) for x in cores]
        for job in cf.as_completed(jobs):
            z=job.result();done.append((z['T'],z['N']));total.update(z['counts']);
            print(json.dumps({'done':len(done),'expected':len(cores),'T':z['T'],'N':z['N'],'row4':z['row4_states']}),flush=True)
            if time.monotonic()>deadline: raise TimeoutError('global time cap')
    status='complete' if len(cores)==1140 else 'partial_scope'
    result={'status':status,'scope':'Necessary S={1}, degree8 rows2--4 only; no d/e/f/g/root screen or arbitrary-cofactor exclusion.','code_sha256':code_hash,'core_list_sha256':core_list_sha256,'cores':len(cores),'completed':len(done),'unknown_cores':len(cores)-len(done),'counts':dict(total),'row4_states':total['row4_states'],'workers':args.workers,'seconds':time.monotonic()-start,'python':sys.version}
    atomic_json(args.out/'summary.json',result);print(json.dumps(result,sort_keys=True))


if __name__=='__main__': main()
