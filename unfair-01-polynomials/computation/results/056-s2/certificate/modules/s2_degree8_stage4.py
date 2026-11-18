#!/usr/bin/env python3
"""Exact necessary rows 1--4 screen for degree-8 quadratic S={2}.

P=1+Lx+a x^2+b x^3+c x^4+d x^5+e x^6+f x^7+x^8,
a^2=T*a-N, with a=(T+sqrt(Delta))/2.  This is a necessary-prefix
enumeration only: it neither searches tail coefficients nor proves the
original uniform-sum problem.
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, hashlib, json, math, os, sys, time
from collections import Counter
from pathlib import Path

Pair=tuple[int,int]                 # p+q*a
ZERO=(0,0); ONE=(1,0)
B2=(17,1); B3=(65,2); NEWTON4_BOUND=(257,4)

def add(x:Pair,y:Pair)->Pair:return (x[0]+y[0],x[1]+y[1])
def scale(k:int,x:Pair)->Pair:return (k*x[0],k*x[1])
def mul(x:Pair,y:Pair,T:int,N:int)->Pair:
    p,q=x; r,s=y
    return (p*r-q*s*N,p*s+q*r+q*s*T)
def conj(x:Pair,T:int)->Pair:
    p,q=x; return (p+q*T,-q)

def sign_plus(x:Pair,T:int,D:int)->int:
    """Exact sign of p+q*(T+sqrt(D))/2; all caller values are nonzero."""
    p,q=x; r=2*p+q*T
    if q==0:return (r>0)-(r<0)
    if q>0 and r>=0:return 1
    if q<0 and r<=0:return -1
    cmp=(r*r>q*q*D)-(r*r<q*q*D)
    # r+q sqrt(D): for q>0,r<0 the sign is opposite sign(r^2-q^2D);
    # for q<0,r>0 it is the same.
    return -cmp if q>0 else cmp

def both_abs(x:Pair, bound:tuple[int,int],T:int,D:int)->bool:
    """Strict |x|<bound in both real embeddings, without floating point."""
    num,den=bound; r=abs(2*x[0]+x[1]*T); q=abs(x[1]); limit=2*num
    gap=limit-r*den
    return gap>0 and q*q*D*den*den < gap*gap

def floor_plus(x:Pair,T:int,D:int)->int:
    p,q=x; radius=abs(p)+abs(q)*(abs(T)+math.isqrt(D)+3)+4
    lo,hi=-radius,radius
    while hi-lo>1:
        mid=(lo+hi)//2
        if sign_plus((p-mid,q),T,D)>=0:lo=mid
        else:hi=mid
    return lo

def unit_shifts(base:Pair,T:int,D:int):
    """All integer h for base+h in [0,1], retaining both rational endpoints."""
    if base[1]==0:
        return [(-base[0],ZERO),(1-base[0],ONE)]
    h=-floor_plus(base,T,D); z=(base[0]+h,base[1])
    assert sign_plus(z,T,D)>0 and sign_plus(add(ONE,scale(-1,z)),T,D)>0
    return [(h,z)]

def binary(x:Pair)->bool:return x[1]==0 and x[0] in (0,1)
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def atomic(path:Path,obj):
    tmp=path.with_suffix(path.suffix+'.partial'); path.parent.mkdir(parents=True,exist_ok=True)
    with tmp.open('w') as f:
        json.dump(obj,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

def core_list():
    out=[]
    for L in range(-9,10):
      for T in range(-82,83):
        x=17-abs(L*L-T); y=83-abs(T)
        if x<=0 or y<=0:continue
        cap=min(x*x,y*y)
        for D in range(1,cap):
          if math.isqrt(D)**2==D or D%4 != (T*T)%4:continue
          N=(T*T-D)//4
          if both_abs((0,1),(83,2),T,D):out.append((L,T,N,D))
    return sorted(out)

def symbolic_check():
    import sympy as sp
    L,a,b,c,d,e,f=sp.symbols('L a b c d e f')
    G1,G2,G3,G4,G5,G6=sp.symbols('G1:7')
    P=[1,L,a,b,c,d,e,f,1]; G=[1,G1,G2,G3,G4,G5,G6]
    A=[sp.expand(sum(P[i]*G[k-i] for i in range(k+1))) for k in range(7)]
    assert sp.expand(A[2]-(a+G2+L*G1))==0
    assert sp.expand(A[3]-(G3+b+L*G2+a*G1))==0
    assert sp.expand(A[4]-(G4+c+L*G3+a*G2+b*G1))==0
    # In the degree-six low coefficient, the only nonrational contribution
    # is a*G4. Thus the separately established high-end endpoint condition
    # A6=B6=0 forces G4=Q4=0 for integer cofactors.
    assert sp.diff(A[6],a)==G4
    # Newton power sums after a^2=T*a-N.
    T,N=sp.symbols('T N')
    s2=L**2-2*a
    s3=-L**3+3*L*a-3*b
    s4=L**4-4*L**2*a+2*(T*a-N)+4*L*b-4*c
    assert sp.expand(s4-(L**4-2*N+4*L*b-4*c+(-4*L**2+2*T)*a))==0

def one(job):
    L,T,N,D,outdir,codehash,maxseconds=job; began=time.monotonic(); ct=Counter(); final=[]
    target=Path(outdir)/'cores'/f'L{L}_T{T}_N{N}.json'
    binding={'version':1,'core':[L,T,N,D],'bounds':{'a_house':[83,2],'s2':[17,1],'s3':[65,2],'s4':[257,4],'b':[-92,92],'c':[-120,120]},'code_sha256':codehash}
    if target.exists():
        x=json.loads(target.read_text())
        if x.get('binding')==binding and x.get('status')=='complete':return x['summary']
        raise RuntimeError(f'refuse stale/nonmatching cache {target}')
    for A1,B1 in ((0,0),(1,0),(0,1)):
      G1=A1-L; Q1=B1-L
      for rho,A2 in unit_shifts((0,1),T,D):
       for sigma,B2 in unit_shifts(conj((0,1),T),T,D):
        c2=add(add(A2,B2),scale(A1*B1,ONE))
        if not binary(c2):ct['c2_reject']+=1;continue
        ct['row2_states']+=1
        G2=rho-L*G1; Q2=sigma-L*Q1
        for b in range(-92,93):
          ct['b_box']+=1
          s3=(-L**3-3*b,3*L)
          if not both_abs(s3,B3,T,D):ct['s3_reject']+=1;continue
          ct['s3_pass']+=1
          U3=(L*G2,G1); V3=add((L*Q2,0),conj((0,Q1),T))
          for eta,A3 in unit_shifts(U3,T,D):
           for zeta,B3v in unit_shifts(V3,T,D):
            c3=add(add(add(A3,B3v),scale(A1,B2)),scale(B1,A2))
            if not binary(c3):ct['c3_reject']+=1;continue
            ct['row3_states']+=1
            G3=eta-b; Q3=zeta-b
            W4=add(add((L*G3+b*G1,G2),ZERO),ZERO)
            Z4=add((L*Q3+b*Q1,0),conj((0,Q2),T))
            # G4=Q4=0 are imposed high-end necessary constraints; c is the
            # same rational coefficient on both sides, so intersect unit shifts.
            left={h:x for h,x in unit_shifts(W4,T,D) if -120<=h<=120}
            right={h:x for h,x in unit_shifts(Z4,T,D) if -120<=h<=120}
            for c in sorted(set(left)&set(right)):
              A4,B4=left[c],right[c]
              s4=(L**4-2*N+4*L*b-4*c,-4*L*L+2*T)
              if not both_abs(s4,NEWTON4_BOUND,T,D):ct['s4_reject']+=1;continue
              ct['s4_pass']+=1
              c4=add(add(add(add(A4,B4),scale(A1,B3v)),mul(A2,B2,T,N)),scale(B1,A3))
              if not binary(c4):ct['c4_reject']+=1;continue
              ct['row4_states']+=1
              final.append([L,T,N,D,b,c,G1,Q1,rho,sigma,eta,zeta])
        if time.monotonic()-began>maxseconds:raise TimeoutError(f'core cap {L,T,N,D}')
    summary={'L':L,'T':T,'N':N,'Delta':D,'counts':dict(ct),'row4_states':len(final),'seconds':time.monotonic()-began}
    atomic(target,{'status':'complete','binding':binding,'summary':summary,'records':final})
    return summary

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--workers',type=int,default=1); ap.add_argument('--max-seconds',type=float,default=14_400); ap.add_argument('--limit',type=int); a=ap.parse_args()
    symbolic_check(); cores=core_list(); codehash=sha(Path(__file__))
    if a.limit is not None: cores=cores[:a.limit]
    jobs=[(*x,str(a.out),codehash,a.max_seconds) for x in cores]
    results=[]
    with cf.ProcessPoolExecutor(max_workers=a.workers) as ex:
      for r in ex.map(one,jobs):results.append(r)
    ct=Counter()
    for r in results:ct.update(r['counts'])
    h=hashlib.sha256()
    for x in cores:h.update(json.dumps(x,separators=(',',':')).encode()+b'\n')
    atomic(a.out/'summary.json',{'status':'complete','scope':'Necessary degree-8 S={2} rows1--4 only; no tail search, root screen, or general conjecture conclusion.','code_sha256':codehash,'cores':len(cores),'completed':len(results),'core_list_sha256':h.hexdigest(),'counts':dict(ct),'row4_states':sum(r['row4_states'] for r in results),'workers':a.workers,'python':sys.version})
    print(json.dumps({'status':'complete','cores':len(cores),'row4_states':sum(r['row4_states'] for r in results),'counts':dict(ct)},sort_keys=True))
if __name__=='__main__':main()
