#!/usr/bin/env python3
"""Independent integral-basis S12 replay, enumerating d, f and g explicitly."""

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

    def integer_box(base,k,radius):
        num,den=BOUND[k]
        m,n=base
        margin=2*num-isqrt(D*(den*radicand_factor*n)**2)-1
        if margin<0:
            return range(0)
        rad=margin//den
        u=2*m+t*n
        return range(max(-radius,(u-rad+2*k-1)//(2*k)),min(radius,(u+rad)//(2*k))+1)

    def desc_norm(P):
        H=[0]*17
        Pb=list(map(conjugate,P))
        for i,a in enumerate(P):
            for j,b in enumerate(Pb):
                z=mul(a,b)
                H[i+j]+=z[0]
        norm_check=[sum(mul(P[i],Pb[k-i])[1] for i in range(max(0,k-8),min(8,k)+1))
                    for k in range(17)]
        assert not any(norm_check)
        return H

    def open_disk_asc(a):
        n=len(a)-1
        while n:
            at,an=a[0],a[n]
            if abs(at)>=abs(an):
                return n
            reduced=[an*a[j]-at*a[n-j] for j in range(1,n+1)]
            f=gcd(*reduced)
            a=[x//f for x in reduced]
            n-=1
        return 0

    def root_code(P):
        for point in ((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6)):
            p,q=point
            val=Z
            for j,x in enumerate(P):
                val=add(val,mulc(p**j*q**(8-j),x))
            for side,x in enumerate((val,conjugate(val))):
                if sign(x)<=0:
                    return "root_positive_reject",{"point":point,"side":side,"value_basis":x}
        h=desc_norm(P)
        failure=open_disk_asc([c*2**i for i,c in enumerate(h)])
        if failure:
            return "root_outer_reject",{"degree":failure,"H":h}
        failure=open_disk_asc([c*2**i for i,c in enumerate(reversed(h))])
        if failure:
            return "root_inner_reject",{"degree":failure,"H":h}
        return "root_pass",{"H":h}

    counts=Counter()
    cases={}
    stage5=[]
    all_a=[(m,n) for n in range(1,15) for m in range(-14,15) if both_abs((m,n),10,1)]
    if not all_a:
        return counts,cases,stage5
    all_b=[(u,w) for w in range(-60,61) if w for u in range(-60,61) if both_abs((u,w),83,2)]
    for a in all_a:
        m,r=a
        ab=conjugate(a)
        A1=-choices(a)[0][0]
        g1=-A1
        q1,b1=choices(ab)[0]
        a1=add(a,integer(g1))
        assert add(a1,b1)==I
        A=(I,a1)
        B=(I,b1)
        s1=mulc(-1,a)
        asq=mul(a,a)
        for b in all_b:
            u,s=b
            bb=conjugate(b)
            counts["raw_core"]+=1
            s2=add(asq,mulc(-2,b))
            if not both_abs(s2,*BOUND[2]):
                continue
            counts["newton2_core"]+=1
            row2=update(A,B,add(b,mulc(g1,a)),add(bb,mulc(q1,ab)),2,(I,))
            counts["row2"]+=len(row2)
            for g2,q2,a2,b2 in row2:
                gs,qs=(1,g1,g2),(1,q1,q2)
                A2,B2=A+(a2,),B+(b2,)
                P2=[I,a,b]
                S2=[integer(8),s1,s2]
                row3=update(A2,B2,engine(P2+[Z],gs,3),engine([I,ab,bb,Z],qs,3),3,(I,))
                counts["row3_normalized"]+=len(row3)
                for g3base,q3base,a3,b3 in row3:
                    A3,B3=A2+(a3,),B2+(b3,)
                    b3root=newton_base(P2,S2,3)
                    for c in integer_box(b3root,3,92):
                        counts["newton3"]+=1
                        P3=P2+[integer(c)]
                        Q3=list(map(conjugate,P3))
                        S3=S2+[add(b3root,integer(-3*c))]
                        G3,Q3coef=gs+(g3base-c,),qs+(q3base-c,)
                        row4=update(A3,B3,engine(P3+[Z],G3,4),engine(Q3+[Z],Q3coef,4),4,(I,))
                        counts["row4_normalized"]+=len(row4)
                        b4root=newton_base(P3,S3,4)
                        for g4base,q4base,a4,b4 in row4:
                            A4,B4=A3+(a4,),B3+(b4,)
                            if (g4base-q4base)*s%r==0:
                                counts["row4_difference_divisibility"]+=1
                            for d in integer_box(b4root,4,120):
                                g4,q4=g4base-d,q4base-d
                                if s*g4%r or s*q4%r: continue
                                g5,q5=-s*g4//r,-s*q4//r
                                G4,Q4=G3+(g4,),Q3coef+(q4,)
                                P4,Q4poly=P3+[integer(d)],Q3+[integer(d)]
                                S4=S3+[add(b4root,integer(-4*d))]
                                v5=add(engine(P4+[Z],G4,5),integer(g5))
                                w5=add(engine(Q4poly+[Z],Q4,5),integer(q5))
                                cross5=cross(A4,B4,5)
                                for e,a5 in choices(v5):
                                    if not -92<=e<=92: continue
                                    b5=add(w5,integer(e))
                                    if sign(b5)<0 or sign(add(I,mulc(-1,b5)))<0: continue
                                    if add(add(a5,b5),cross5) not in (Z,I): continue
                                    s5=add(newton_base(P4,S4,5),integer(-5*e))
                                    if not both_abs(s5,*BOUND[5]): continue
                                    counts["row5"]+=1
                                    A5,B5=A4+(a5,),B4+(b5,)
                                    if cross(A5,B5,6) not in (Z,I): continue
                                    counts["row6_cross_pass"]+=1
                                    state=(D,2*m+t*r,r,2*u+t*s,s,c,d,e)
                                    stage5.append(state)
                                    A6,B6=A5+(Z,),B5+(Z,)
                                    G5,Q5=G4+(g5,),Q4+(q5,)
                                    P5,PB5=P4+[integer(e)],Q4poly+[integer(e)]
                                    S5=S4+[s5]
                                    b6root=newton_base(P5,S5,6)
                                    for f in integer_box(b6root,6,41):
                                        P6,PB6=P5+[integer(f)],PB5+[integer(f)]
                                        s6=add(b6root,integer(-6*f))
                                        S6=S5+[s6]
                                        aa6,bb6=engine(P6,G5,6),engine(PB6,Q5,6)
                                        assert aa6[1]==bb6[1]==0
                                        G6,Q6=G5+(-aa6[0],),Q5+(-bb6[0],)
                                        b7root=newton_base(P6,S6,7)
                                        for g in integer_box(b7root,7,9):
                                            P7,PB7=P6+[integer(g)],PB6+[integer(g)]
                                            s7=add(b7root,integer(-7*g))
                                            s8=add(newton_base(P7,S6+[s7],8),integer(-8))
                                            if not both_abs(s8,*BOUND[8]): continue
                                            baseA,baseB=engine(P7,G6,7),engine(PB7,Q6,7)
                                            row7=update(A6,B6,baseA,baseB,7)
                                            for g7,q7,a7,b7 in row7:
                                                if r!=s*g and (a7!=Z or b7!=Z): continue
                                                counts["row7"]+=1
                                                A7,B7=A6+(a7,),B6+(b7,)
                                                P8,PB8=P7+[I],PB7+[I]
                                                r8=update(A7,B7,engine(P8,G6+(g7,),8),
                                                          engine(PB8,Q6+(q7,),8),8)
                                                if r8:
                                                    counts["row8_states"]+=len(r8)
                                                    row=(*state,f,g)
                                                    if row not in cases:
                                                        cases[row]=root_code(P8)
    counts["mothers"]=len(cases)
    counts.update(v[0] for v in cases.values())
    return counts,cases,stage5


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    a=ap.parse_args()
    assert not a.out.exists()
    domains=[D for D in range(2,400) if all(D%(p*p) for p in range(2,isqrt(D)+1))]
    totals=Counter()
    rows={}
    states=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for counts,local,stage5 in pool.map(audit_domain,domains,chunksize=1):
            totals.update(counts)
            assert rows.keys().isdisjoint(local)
            rows.update(local)
            states.extend(stage5)
    imported=a.reference.read_bytes()
    expected={tuple(row):code for row,code in map(json.loads,imported.splitlines())}
    assert len(expected)==len(imported.splitlines())
    assert {row:code for row,(code,witness) in rows.items()}==expected
    reference_totals=json.loads(a.reference.with_name("full-entry.json").read_bytes())["counts"]
    assert all(value==reference_totals[key] for key,value in totals.items())
    out=dict(status="verified",counts=dict(totals),reference_sha256=sha256(imported).hexdigest(),
             stage5=sorted(states),rows=[(row,*data) for row,data in sorted(rows.items())])
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(dict(status="verified",counts=dict(totals),stage5=sorted(states)),sort_keys=True))


if __name__=="__main__":
    main()
