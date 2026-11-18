#!/usr/bin/env python3
"""Independent integral-basis S123 replay using complete integral-basis coordinates and direct common shifts."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path

Z=(0,0)
I=(1,0)
BOUND={2:(17,1),3:(65,2),4:(257,4),5:(1025,8),6:(4097,16),7:(16385,32),8:(65537,64)}


def add(x,y): return x[0]+y[0],x[1]+y[1]
def mulc(k,x): return k*x[0],k*x[1]
def integer(x): return x,0


def audit_domain(D):
    t=int(D%4==1)
    radicand_factor=1 if t else 2
    omega_constant=(D-1)//4 if t else D

    def mul(x,y):
        p,q=x
        r,s=y
        return p*r+q*s*omega_constant,p*s+q*r+q*s*t

    def conjugate(x):
        return x[0]+t*x[1],-x[1]

    def sign(x):
        m,n=x
        u=2*m+t*n
        w=radicand_factor*n
        if u>=0 and w>=0: return int(u+w>0)
        if u<=0 and w<=0: return -int(u+w<0)
        q=u*u-D*w*w
        assert q
        return int(q>0)*2-1 if u>0 else 1-int(q>0)*2

    def both_abs(x,num,den):
        m,n=x
        gap=2*num-den*abs(2*m+t*n)
        return gap>0 and gap*gap>D*(den*radicand_factor*n)**2

    def choices(x):
        m,n=x
        if not n:
            return ((-m,Z),(1-m,I))
        R=abs(m)+abs(n)*(D+1)+2
        lo,hi=-R,R
        while hi-lo>1:
            mid=(lo+hi)//2
            if sign((m-mid,n))>=0: lo=mid
            else: hi=mid
        return ((-lo,(m-lo,n)),)

    def engine(P,G,k):
        out=Z
        for i in range(1,k+1):
            out=add(out,mulc(G[k-i],P[i]))
        return out

    def cross(A,B,k):
        out=Z
        for i in range(1,k):
            out=add(out,mul(A[i],B[k-i]))
        return out

    def update(A,B,u,w,k,req=(Z,I)):
        h=cross(A,B,k)
        return [(g,q,a,b) for g,a in choices(u) for q,b in choices(w)
                if add(add(a,b),h) in req]

    def newton_base(P,S,k):
        base=Z
        for j in range(1,k):
            base=add(base,mulc(-1,mul(P[j],S[k-j])))
        return base

    totals=Counter()
    witnesses=[]
    all_b=[(m,s) for s in range(-38,39) if s for m in range(-61,62) if both_abs((m,s),83,2)]
    for r in range(1,10):
        for m in range(-14,15):
            a=(m,r)
            if not both_abs(a,10,1):continue
            totals["a_cores"]+=1
            ab=conjugate(a)
            (g1,a1),=choices(a)
            (q1,b1),=choices(ab)
            assert add(a1,b1)==I
            A,B=(I,a1),(I,b1)
            aa=mul(a,a)
            for b in all_b:
                u,s=b
                bb=conjugate(b)
                S2=add(aa,mulc(-2,b))
                if not both_abs(S2,17,1):continue
                P2=[I,a,b]
                S=[integer(8),mulc(-1,a),S2]
                for g in range(-9,10):
                    if not g or s%g:continue
                    t_=s//g
                    totals["newton2_ab"]+=1
                    for g2,q2,a2,b2 in update(A,B,add(b,mulc(g1,a)),add(bb,mulc(q1,ab)),2,(I,)):
                        totals["row2"]+=1
                        A2,B2=A+(a2,),B+(b2,)
                        c0,cb0=(0,t_),conjugate((0,t_))
                        U3=add(c0,add(mulc(g1,b),mulc(g2,a)))
                        V3=add(cb0,add(mulc(q1,bb),mulc(q2,ab)))
                        for rho,sigma,a3,b3 in update(A2,B2,U3,V3,3,(I,)):
                            totals["row3_normalized"]+=1
                            numerator=s*(rho-sigma)+t_*(g2-q2)
                            if numerator%r:continue
                            totals["row5_difference_divisibility"]+=1
                            base3=add(newton_base(P2,S,3),mulc(-3,c0))
                            trace_base=2*base3[0]+t*base3[1]
                            vlo=-((64-trace_base)//6)
                            vhi=(64+trace_base)//6
                            for v in range(vlo,vhi+1):
                                c=(v,t_)
                                if not both_abs(c,185,2):continue
                                S3=add(base3,integer(-3*v))
                                if not both_abs(S3,65,2):continue
                                cb=conjugate(c)
                                G,Q=(1,g1,g2,rho-v),(1,q1,q2,sigma-v)
                                A3,B3=A2+(a3,),B2+(b3,)
                                P,Pb=P2+[c,Z],[I,ab,bb,cb,Z]
                                baseA,baseB=engine(P,G,4),engine(Pb,Q,4)
                                diff=add(add(baseA,mulc(-1,baseB)),integer(-numerator//r))
                                if sign(add(I,diff))<0 or sign(add(I,mulc(-1,diff)))<0:
                                    continue
                                totals["row4_difference"]+=1
                                left,right=s*G[3]+t_*g2,s*Q[3]+t_*q2
                                if left%r or right%r:continue
                                totals["row5_full_divisibility"]+=1
                                U4,V4=add(baseA,integer(-left//r)),add(baseB,integer(-right//r))
                                cross4=cross(A3,B3,4)
                                val=add(add(U4,V4),cross4)
                                assert val[1]==0
                                if val[0]%2:
                                    totals["row4_parity"]+=1
                                for d,a4 in choices(U4):
                                    b4=add(V4,integer(d))
                                    if sign(b4)<0 or sign(add(I,mulc(-1,b4)))<0:
                                        continue
                                    totals["shared_unit_d"]+=1
                                    C4=add(add(a4,b4),cross4)
                                    assert C4[1]==0
                                    totals["C4_value_"+str(C4[0])]+=1
                                    row=(D,2*m+t*r,r,2*u+t*s,s,2*v+t*t_,t_,d,g)
                                    witnesses.append((row,C4[0],a4,b4))
                                    if C4 in (Z,I):
                                        totals["shared_d_binary"]+=1
    return totals,witnesses


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,400) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    total=Counter()
    rows=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for ct,rs in pool.map(audit_domain,domains,chunksize=1):
            total.update(ct)
            rows.extend(rs)
    reference=a.reference.read_bytes()
    prior=json.loads(reference)["counts"]
    for k in ("a_cores","newton2_ab","row2","row3_normalized","row5_difference_divisibility",
              "row4_difference","row5_full_divisibility","row4_parity","shared_d_binary"):
        assert total[k]==prior.get(k,0), (k,total[k],prior.get(k,0))
    assert total["shared_d_binary"]==0
    rows=sorted(rows)
    payload=b"".join(json.dumps(row,separators=(",",":")).encode()+b"\n" for row in rows)
    a.out.with_suffix(".rows.jsonl").write_bytes(payload)
    out=dict(status="verified",counts=dict(total),rows_sha256=sha256(payload).hexdigest(),
             reference_sha256=sha256(reference).hexdigest())
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True))


if __name__=="__main__":
    main()
