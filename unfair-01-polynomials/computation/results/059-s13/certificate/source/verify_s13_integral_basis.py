#!/usr/bin/env python3
"""Independent integral-basis S13 replay using complete integral-basis coordinates and direct common shifts."""

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
    v=(D-1)//4 if t else D

    def mul(x,y):
        p,q=x
        r,s=y
        return p*r+q*s*v,p*s+q*r+q*s*t

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

    counts=Counter()
    contradictions=[]
    a_values=[(m,r) for r in range(1,15) for m in range(-14,15) if both_abs((m,r),10,1)]
    if not a_values:
        return counts,contradictions
    b_values=[(u,s) for s in range(-84,85) if s for u in range(-134,135)
              if both_abs((u,s),185,2)]
    for a in a_values:
        m,r=a
        ab=conjugate(a)
        (g1,a1),=(choices(a))
        (q1,b1),=(choices(ab))
        assert add(a1,b1)==I
        A,B=(I,a1),(I,b1)
        counts["a_cores"]+=1
        row2=update(A,B,mulc(g1,a),mulc(q1,ab),2,(I,))
        assert len(row2)<=1
        if not row2:
            continue
        counts["a_row2"]+=1
        rho2,sigma2,a2,b2=row2[0]
        A2,B2=A+(a2,),B+(b2,)
        asq=mul(a,a)
        atimes=mul(a,asq)
        for c in range(-41,42):
            s2=add(asq,integer(-2*c))
            if not both_abs(s2,*BOUND[2]):
                continue
            counts["newton2_c"]+=1
            g2,q2=rho2-c,sigma2-c
            k3=add(mulc(-1,atimes),mulc(3*c,a))
            for b in b_values:
                u,s=b
                s3=add(k3,mulc(-3,b))
                if not both_abs(s3,*BOUND[3]):
                    continue
                counts["newton3_core"]+=1
                if s*g2%r or s*q2%r:
                    continue
                counts["row5_projection_divisibility"]+=1
                bb=conjugate(b)
                left=add(add(b,mulc(g1,integer(c))),mulc(g2,a))
                right=add(add(bb,mulc(q1,integer(c))),mulc(q2,ab))
                row3=update(A2,B2,left,right,3,(I,))
                assert len(row3)<=1
                if not row3:
                    continue
                counts["row3"]+=1
                g3,q3,a3,b3=row3[0]
                A3,B3=A2+(a3,),B2+(b3,)
                g4,q4=-s*g2//r,-s*q2//r
                p4_base,q4_base=engine([I,a,integer(c),b,Z],(1,g1,g2,g3),4),engine(
                    [I,ab,integer(c),bb,Z],(1,q1,q2,q3),4)
                row4=update(A3,B3,p4_base,q4_base,4,(I,))
                counts["row4_normalized"]+=len(row4)
                shift_a,shift_b=add(p4_base,integer(g4)),add(q4_base,integer(q4))
                for d,a4 in choices(shift_a):
                    b4=add(shift_b,integer(d))
                    if sign(b4)<0 or sign(add(I,mulc(-1,b4)))<0:
                        continue
                    counts["shared_unit_d"]+=1
                    case=(D,2*m+t*r,r,c,2*u+t*s,s,d)
                    total4=add(add(a4,b4),cross(A3,B3,4))
                    assert total4[1]==0
                    contradictions.append((case,total4[0],a4,b4))
                    counts["C4_value_"+str(total4[0])]+=1
                    if total4 in (Z,I):
                        counts["shared_d"]+=1
                        if -120<=d<=120:
                            counts["shared_d_box"]+=1
                            s4=add(newton_base([I,a,integer(c),b],[integer(8),mulc(-1,a),s2,s3],4),
                                   integer(-4*d))
                            if both_abs(s4,*BOUND[4]):
                                counts["row4"]+=1
    return counts,contradictions


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,400) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    totals=Counter()
    rows=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for counts,local in pool.map(audit_domain,domains,chunksize=1):
            totals.update(counts)
            rows.extend(local)
    expected=json.loads(a.reference.read_bytes())
    for k in ("a_cores","a_row2","newton2_c","newton3_core","row5_projection_divisibility",
              "row3","row4_normalized","shared_d","shared_d_box","row4"):
        assert totals[k]==expected["counts"].get(k,0), (k,totals[k],expected["counts"])
    assert totals["shared_d"]==totals["row4"]==0
    rows.sort()
    payload=b"".join(json.dumps(r,separators=(",",":")).encode()+b"\n" for r in rows)
    a.out.with_suffix(".rows.jsonl").write_bytes(payload)
    result=dict(status="verified",counts=dict(totals),reference_sha256=sha256(a.reference.read_bytes()).hexdigest(),
                rows_sha256=sha256(payload).hexdigest())
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
