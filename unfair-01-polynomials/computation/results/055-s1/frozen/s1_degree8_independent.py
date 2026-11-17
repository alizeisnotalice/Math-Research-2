#!/Users/zzh/.venvs/main/bin/python
"""Independent exact necessary-prefix audit for degree-8 quadratic S={1}.

This implementation deliberately uses pairs (r,s) for (r+s*sqrt(D))/2,
not the producing program's p+q*a representation.  It uses direct A*B
convolutions at every retained row and never imports the producing programs.

Scope: complete finite necessary boxes through rows 2--7, including the
high-end constraints G6=Q6=0 and A7=B7=0.  Surviving prefixes are not
factorizations or counterexamples.
"""
from __future__ import annotations

import argparse
import concurrent.futures as futures
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

# (r,s) denotes (r+s sqrt(D))/2.  Integral n is (2n,0).
Pair = tuple[int, int]
Z: Pair = (0, 0)
ONE: Pair = (2, 0)
BOUNDS = {2:(17,1), 3:(65,2), 4:(257,4), 5:(1025,8),
          6:(4097,16), 7:(16385,32), 8:(65537,64)}


def I(n: int) -> Pair: return (2*n, 0)
def add(x: Pair, y: Pair) -> Pair: return (x[0]+y[0], x[1]+y[1])
def neg(x: Pair) -> Pair: return (-x[0], -x[1])
def scale(n: int, x: Pair) -> Pair: return (n*x[0], n*x[1])
def conj(x: Pair) -> Pair: return (x[0], -x[1])


def mul(x: Pair, y: Pair, D: int) -> Pair:
    r,s=x; u,v=y
    a,b=r*u+s*v*D, r*v+s*u
    assert a % 2 == 0 and b % 2 == 0, (x,y,D,a,b)
    return (a//2,b//2)


def sign(x: Pair, D: int) -> int:
    """Sign of (r+s sqrt(D))/2, using only integer squares."""
    r,s=x
    if s == 0: return (r > 0) - (r < 0)
    if r >= 0 and s > 0: return 1
    if r <= 0 and s < 0: return -1
    q = s*s*D
    p = r*r
    if r > 0:  # s<0
        return (p > q) - (p < q)
    # r<0 and s>0
    return (q > p) - (q < p)


def floor_q(x: Pair, D: int) -> int:
    r,s=x
    bound = abs(r) + abs(s)*(math.isqrt(D)+2) + 4
    lo,hi=-bound,bound
    while hi-lo > 1:
        m=(lo+hi)//2
        if sign(add(x,I(-m)),D) >= 0: lo=m
        else: hi=m
    return lo


def unit_shifts(base: Pair, D: int) -> list[tuple[int, Pair]]:
    """All integral h with base+h in [0,1]."""
    if base[1] == 0:
        assert base[0] % 2 == 0
        n=base[0]//2
        return [(-n,Z), (1-n,ONE)]
    h=-floor_q(base,D); out=add(base,I(h))
    assert sign(out,D)>0 and sign(add(ONE,neg(out)),D)>0
    return [(h,out)]


def both_abs_lt(x: Pair, num: int, den: int, D: int) -> bool:
    # |x|<num/den under both embeddings.
    for y in (x,conj(x)):
        if sign(add(I(num),neg(scale(den,y))),D) <= 0: return False
        if sign(add(I(num),scale(den,y)),D) <= 0: return False
    return True


def coeff(P: list[Pair], G: list[Pair], k: int, D: int) -> Pair:
    z=Z
    for i in range(k+1):
        if i < len(P) and k-i < len(G): z=add(z,mul(P[i],G[k-i],D))
    return z


def cbit(A: list[Pair], B: list[Pair], k: int, D: int) -> int | None:
    z=Z
    for i in range(k+1): z=add(z,mul(A[i],B[k-i],D))
    if z == Z: return 0
    if z == ONE: return 1
    return None


def newton(P: list[Pair], upto: int, D: int) -> list[Pair]:
    """For P=prod(1-w_i x), s_k=-sum_{i<k}p_i*s_{k-i}-k p_k."""
    s=[Z]
    for k in range(1,upto+1):
        z=Z
        for i in range(1,k): z=add(z,mul(P[i],s[k-i],D))
        s.append(neg(add(z,scale(k,P[k]))))
    return s


def cores() -> list[tuple[int,int,int]]:
    ans=[]
    for T in range(-19,20):
        for N in range(-99,100):
            D=T*T-4*N
            if D<=0 or math.isqrt(D)**2==D: continue
            if D >= (20-abs(T))**2: continue
            ans.append((T,N,D))
    assert len(ans)==1140
    return ans


def atom(path: Path, obj: object) -> None:
    q=path.with_suffix(path.suffix+'.partial')
    with q.open('w') as f:
        json.dump(obj,f,sort_keys=True,separators=(',',':')); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(q,path)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage_core(core: tuple[int,int,int], outdir: str, code_sha: str, seconds: float) -> dict:
    T,N,D=core; target=Path(outdir)/'cores'/f'T{T}_N{N}.json'
    binding={'code_sha256':code_sha,'core':[T,N,D], 'representation':'(r+s*sqrt(D))/2',
             'boxes':{'b':[-41,41],'c':[-92,92],'d':[-120,120],'e':[-92,92],'f':[-41,41],'g':[-9,9]}}
    if target.exists():
        old=json.loads(target.read_text()); assert old['binding']==binding and old['status']=='complete'
        return old['summary'] | {'cached':True}
    began=time.monotonic(); a=(T,1); abar=(T,-1)
    g1=-floor_q(a,D); q1=-floor_q(abar,D)
    A1=add(a,I(g1)); B1=add(abar,I(q1))
    assert sign(A1,D)>0 and sign(add(ONE,neg(A1)),D)>0
    assert sign(B1,D)>0 and sign(add(ONE,neg(B1)),D)>0 and cbit([ONE,A1],[ONE,B1],1,D)==1
    ct=Counter(cores=1,b_box=83); row4=[]
    # The loops retain every current direct-AB prefix, then extend it.
    for b in range(-41,42):
        P2=[ONE,a,I(b)]
        s2=newton(P2,2,D)[2]
        if not both_abs_lt(s2,*BOUNDS[2],D): ct['s2_reject']+=1; continue
        ct['s2_pass']+=1
        # Use rho=G2+b and sigma=Q2+b as the row-2 integral shifts.
        baseA2=scale(g1,a); baseB2=scale(q1,abar)
        for rho,A2 in unit_shifts(baseA2,D):
            for sigma,B2 in unit_shifts(baseB2,D):
                if cbit([ONE,A1,A2],[ONE,B1,B2],2,D) is None: ct['c2_reject']+=1; continue
                ct['row2_states']+=1
                G2=rho-b; Q2=sigma-b
                for c in range(-92,93):
                    P3=[ONE,a,I(b),I(c)]
                    s3=newton(P3,3,D)[3]
                    if not both_abs_lt(s3,*BOUNDS[3],D): ct['s3_reject']+=1; continue
                    ct['s3_pass']+=1
                    GA=[ONE,I(g1),I(G2)]; QB=[ONE,I(q1),I(Q2)]
                    # eta=G3+c and zeta=Q3+c are the row-3 shifts.
                    baseA3=add(coeff(P3,GA,3,D),I(-c)); baseB3=add(coeff([ONE,abar,I(b),I(c)],QB,3,D),I(-c))
                    for eta,A3 in unit_shifts(baseA3,D):
                        for zeta,B3 in unit_shifts(baseB3,D):
                            if cbit([ONE,A1,A2,A3],[ONE,B1,B2,B3],3,D) is None: ct['c3_reject']+=1; continue
                            ct['row3_states']+=1
                            G3=eta-c; Q3=zeta-c
                            for d in range(-120,121):
                                P4=[ONE,a,I(b),I(c),I(d)]
                                s4=newton(P4,4,D)[4]
                                if not both_abs_lt(s4,*BOUNDS[4],D): ct['s4_reject']+=1; continue
                                ct['s4_pass']+=1
                                GA4=[ONE,I(g1),I(G2),I(G3)]; QB4=[ONE,I(q1),I(Q2),I(Q3)]
                                # alpha=G4+d+b*rho-b^2 (and beta) are row-4 shifts.
                                shiftA4=d+b*rho-b*b; shiftB4=d+b*sigma-b*b
                                baseA4=add(coeff(P4,GA4,4,D),I(-shiftA4)); baseB4=add(coeff([ONE,abar,I(b),I(c),I(d)],QB4,4,D),I(-shiftB4))
                                for alpha,A4 in unit_shifts(baseA4,D):
                                    for beta,B4 in unit_shifts(baseB4,D):
                                        if cbit([ONE,A1,A2,A3,A4],[ONE,B1,B2,B3,B4],4,D) is None: ct['c4_reject']+=1; continue
                                        ct['row4_states']+=1
                                        G4=alpha-shiftA4; Q4=beta-shiftB4
                                        # Normalized to the producing schema solely for later comparison.
                                        row4.append([b,c,d,rho,sigma,eta,zeta,alpha,beta])
        if time.monotonic()-began > seconds: raise TimeoutError(f'core cap {core}')

    final=[]
    # Rebuild every retained row4 prefix as direct polynomials, then rows5--7.
    for b,c,d,rho,sigma,eta,zeta,alpha,beta in row4:
        G2=rho-b; Q2=sigma-b; G3=eta-c; Q3=zeta-c; G4=alpha-d-b*rho+b*b; Q4=beta-d-b*sigma+b*b
        P4=[ONE,a,I(b),I(c),I(d)]; Pb4=[ONE,abar,I(b),I(c),I(d)]
        G4v=[ONE,I(g1),I(G2),I(G3),I(G4)]; Q4v=[ONE,I(q1),I(Q2),I(Q3),I(Q4)]
        A=[coeff(P4,G4v,k,D) for k in range(5)]; B=[coeff(Pb4,Q4v,k,D) for k in range(5)]
        assert all(cbit(A,B,k,D) is not None for k in range(5))
        baseA5=coeff(P4,G4v,5,D); baseB5=coeff(Pb4,Q4v,5,D)
        for eta5,A5 in unit_shifts(baseA5,D):
            for zeta5,B5 in unit_shifts(baseB5,D):
                AA=A+[A5]; BB=B+[B5]
                if cbit(AA,BB,5,D) is None: ct['c5_reject']+=1; continue
                ct['row5_states']+=1
                for e in range(-92,93):
                    ct['e_box']+=1
                    P5=P4+[I(e)]; s5=newton(P5,5,D)[5]
                    if not both_abs_lt(s5,*BOUNDS[5],D): ct['s5_reject']+=1; continue
                    ct['s5_pass']+=1
                    G5=eta5-e; Q5=zeta5-e
                    G5v=G4v+[I(G5)]; Q5v=Q4v+[I(Q5)]
                    # G6=Q6=0: f is the common polynomial coefficient and sole row-6 integral shift.
                    baseA6=coeff(P5,G5v,6,D); baseB6=coeff([ONE,abar,I(b),I(c),I(d),I(e)],Q5v,6,D)
                    fA={h for h,_ in unit_shifts(baseA6,D)}; fB={h for h,_ in unit_shifts(baseB6,D)}
                    for f in sorted(fA & fB & set(range(-41,42))):
                        P6=P5+[I(f)]; Pb6=[ONE,abar,I(b),I(c),I(d),I(e),I(f)]
                        A6=add(baseA6,I(f)); B6=add(baseB6,I(f))
                        A6v=AA+[A6]; B6v=BB+[B6]
                        s6=newton(P6,6,D)[6]
                        if not both_abs_lt(s6,*BOUNDS[6],D): ct['s6_reject']+=1; continue
                        ct['s6_pass']+=1
                        if cbit(A6v,B6v,6,D) is None: ct['c6_reject']+=1; continue
                        ct['row6_states']+=1
                        cross=Z
                        for i in range(1,7): cross=add(cross,mul(A6v[i],B6v[7-i],D))
                        if cross not in (Z,ONE): ct['c7_reject']+=1; continue
                        ct['cross7_pass']+=1
                        for g in range(-9,10):
                            ct['g_box']+=1
                            P7=P6+[I(g)]; s7=newton(P7,7,D)[7]; s8=newton(P7+[ONE],8,D)[8]
                            if not both_abs_lt(s7,*BOUNDS[7],D): ct['s7_reject']+=1; continue
                            if not both_abs_lt(s8,*BOUNDS[8],D): ct['s8_reject']+=1; continue
                            G7=-g-(b*G5+c*G4+d*G3+e*G2+f*g1)
                            Q7=-g-(b*Q5+c*Q4+d*Q3+e*Q2+f*q1)
                            GA=G5v+[Z,I(G7)]; QB=Q5v+[Z,I(Q7)]
                            assert coeff(P7,GA,7,D)==Z and coeff([ONE,abar,I(b),I(c),I(d),I(e),I(f),I(g)],QB,7,D)==Z
                            ct['row7_states']+=1
                            final.append([b,c,d,e,f,g,rho,sigma,eta,zeta,alpha,beta,eta5,zeta5,G7,Q7])
        if time.monotonic()-began > seconds: raise TimeoutError(f'core cap {core}')
    assert len(row4)==ct['row4_states'] and len(final)==ct['row7_states']
    summary={'T':T,'N':N,'Delta':D,'counts':dict(ct),'row4_states':len(row4),'row7_states':len(final),'seconds':time.monotonic()-began}
    atom(target,{'status':'complete','binding':binding,'summary':summary,'row4_records':row4,'row7_records':final})
    return summary


def symbolic_check() -> None:
    import sympy as sp
    x=sp.symbols('x'); a,b,c,d,e,f,g=sp.symbols('a b c d e f g'); G1,G2,G3,G4,G5,G6,G7=sp.symbols('G1:8')
    P=1+a*x+b*x**2+c*x**3+d*x**4+e*x**5+f*x**6+g*x**7+x**8
    G=1+G1*x+G2*x**2+G3*x**3+G4*x**4+G5*x**5+G6*x**6+G7*x**7
    A=sp.expand(P*G)
    assert sp.expand(A.coeff(x,6)-(a*G5+b*G4+c*G3+d*G2+e*G1+f+G6))==0
    assert sp.expand(A.coeff(x,7)-(a*G6+b*G5+c*G4+d*G3+e*G2+f*G1+g+G7))==0
    # Check the displayed high-end implication's two exact forced zero equations.
    assert sp.expand(A.coeff(x,7).subs(G6,0)- (b*G5+c*G4+d*G3+e*G2+f*G1+g+G7))==0
    print(json.dumps({'status':'symbolic_passed_direct_A_times_P_rows6_7'}))


def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path); ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--seconds',type=float,default=14400); ap.add_argument('--core-seconds',type=float,default=600)
    ap.add_argument('--core-limit',type=int); ap.add_argument('--only-core',nargs=2,type=int); ap.add_argument('--symbolic-check',action='store_true'); args=ap.parse_args()
    if args.symbolic_check: symbolic_check(); return
    if args.out is None: ap.error('--out is required unless --symbolic-check is used')
    if not __debug__: raise RuntimeError('refuse -O')
    assert 1<=args.workers<=4
    args.out.mkdir(parents=True,exist_ok=True); (args.out/'cores').mkdir(exist_ok=True)
    allcores=cores()
    if args.only_core:
        t,n=args.only_core; work=[z for z in allcores if z[:2]==(t,n)]
        assert len(work)==1
    else:
        work=allcores[:args.core_limit] if args.core_limit else allcores
    start=time.monotonic(); deadline=start+args.seconds; code=sha(Path(__file__)); ct=Counter(); done=[]
    with futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        fs=[pool.submit(stage_core,c,str(args.out),code,args.core_seconds) for c in work]
        for job in futures.as_completed(fs):
            z=job.result(); ct.update(z['counts']); done.append((z['T'],z['N'])); print(json.dumps({'done':len(done),'expected':len(work),'core':[z['T'],z['N']],'row4':z['row4_states'],'row7':z['row7_states']}),flush=True)
            if time.monotonic()>deadline: raise TimeoutError('global time cap')
    status='complete' if len(work)==len(allcores) else 'partial_scope'
    out={'status':status,'scope':'Independent degree-8 S={1} necessary finite prefix screen rows2--7 only; no root screen, no arbitrary-cofactor exclusion.','code_sha256':code,'cores':len(work),'completed':len(done),'unknown_cores':len(work)-len(done),'counts':dict(ct),'row4_states':ct['row4_states'],'row7_states':ct['row7_states'],'workers':args.workers,'seconds':time.monotonic()-start,'python':sys.version}
    atom(args.out/'summary.json',out); print(json.dumps(out,sort_keys=True))


if __name__=='__main__': main()
