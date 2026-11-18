#!/usr/bin/env python3
"""Independent degree-nine singleton S3 enumeration and tree using the (1,a) order."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path

Z=(0,0)
I=(1,0)
A=(0,1)
COEFF=(1,10,51,133,212,212,133,51,10,1)
LIM={k:(4*(2**(2*k)+1)+2**k,2**k) for k in range(1,10)}


def add(x,y):return (x[0]+y[0],x[1]+y[1])
def sub(x,y):return (x[0]-y[0],x[1]-y[1])
def scale(c,x):return (c*x[0],c*x[1])
def z(c):return (c,0)


def env(T,D):
    N=(T*T-D)//4

    def conj(x):
        m,n=x
        return m+T*n,-n

    def key(x):
        m,n=x
        return 2*m+T*n,n

    def mul(x,y):
        m,n=x
        k,l=y
        return m*k-N*n*l,m*l+n*k+T*n*l

    def sign(x):
        u,v=key(x)
        if not v:return (u>0)-(u<0)
        if u==0 or (u>0)==(v>0):return (v>0)-(v<0)
        s=u*u-D*v*v
        return ((s>0)-(s<0))*((u>0)-(u<0))

    def house(x,k):
        u,v=key(x)
        h,l=LIM[k]
        margin=2*h-l*abs(u)
        return margin>0 and margin*margin>D*(l*v)**2

    def shifts(x):
        m,n=x
        if not n:return ((-m,Z),(1-m,I))
        bound=abs(m)+abs(n)*(abs(T)+D+1)+1
        left,right=-bound,bound
        while right-left>1:
            centre=(left+right)//2
            if sign((m-centre,n))<0:right=centre
            else:left=centre
        return ((-left,(m-left,n)),)

    def base(p,U,k):
        out=Z
        for j in range(1,min(k,len(p)-1)+1):
            if 0<=k-j<len(U):
                out=add(out,scale(U[k-j],p[j]))
        return out

    def cross(X,Y,k):
        out=Z
        for j in range(1,k):
            out=add(out,mul(X[j],Y[k-j]))
        return out

    def moment(p,S,k):
        value=Z
        for j in range(1,k):
            value=sub(value,mul(p[j],S[k-j]))
        return value

    def stage(p,states,k,force_g=None,force_coeff=None,against=None):
        pb=tuple(map(conj,p))
        outputs=[]
        for G,Q,X,Y in states:
            u,v=base(p,G,k),base(pb,Q,k)
            cc=cross(X,Y,k)
            for g,x in shifts(u):
                if (force_g is not None and g!=force_g) or (force_coeff is not None and x!=force_coeff):
                    continue
                if against is not None and x!=Z and against[1]!=Z:continue
                for c in (Z,I):
                    y=sub(sub(c,cc),x)
                    if (force_coeff is not None and y!=force_coeff) or sign(y)<0 or sign(sub(I,y))<0:
                        continue
                    if against is not None and y!=Z and against[0]!=Z:continue
                    q,s=sub(y,v)
                    if s or (force_g is not None and q!=force_g):continue
                    outputs.append((G+(g,),Q+(q,),X+(x,),Y+(y,)))
        return outputs

    def bounds(base,k):
        u,v=key(base)
        high,scale_=LIM[k]
        radius=(2*high-isqrt((scale_*v)**2*D)-1)//scale_
        if radius<0:return ()
        return range(max(-COEFF[k],-((radius-u)//(2*k))),min(COEFF[k],(radius+u)//(2*k))+1)

    return locals()


def core_rows():
    initial=((1,),(1,),(I,),(I,))
    for L in range(-10,11):
        for M in range(-51,52):
            if abs(L*L-2*M)>17:continue
            base=-L**3+3*L*M
            for T in range((2*base-66+2)//3,(2*base+66)//3+1):
                allowance=67-abs(2*base-3*T)
                for D in range(T*T%4,(allowance*allowance-1)//9+1,4):
                    if not D or isqrt(D)**2==D:continue
                    if abs(T)+isqrt(D)>=268:continue
                    E=env(T,D)
                    p=(I,z(L),z(M),A)
                    pre=[initial]
                    for k in range(1,4):
                        pre=E["stage"](p,pre,k,force_g=0 if k==3 else None)
                    if pre:yield ((L,M,T,D),len(pre))


def ordinary_moments(p,E,check_all=True):
    S=[z(9)]
    for k in range(1,len(p)):
        value=add(E["moment"](p,S,k),scale(-k,p[k]))
        if check_all and not E["house"](value,k):return None
        S.append(value)
    return S


def root_code(p,E):
    for u,v in ((1,1),(1,2),(2,1),(1,4),(4,1),(3,4),(4,3),(6,5),(5,6)):
        x=Z
        for k,c in enumerate(p):
            x=add(x,scale(u**k*v**(9-k),c))
        if E["sign"](x)<=0 or E["sign"](E["conj"](x))<=0:return "root_positive_reject"
    rational=[x[0] for x in p]
    H=[0]*19
    for i,c in enumerate(rational):
        for j,d in enumerate(rational):
            H[i+j]+=c*d
        H[i+3]+=E["T"]*c
    H[6]+=E["N"]
    for name,s in (("outer",H),("inner",H[::-1])):
        b=[c*2**i for i,c in enumerate(s)]
        while len(b)>1:
            a_n,a_0=b[-1],b[0]
            if abs(a_0)>=abs(a_n):return "root_"+name+"_reject"
            reduced=[a_n*x-a_0*y for x,y in zip(b[1:],b[-2::-1])]
            common=gcd(*reduced)
            b=[x//common for x in reduced]
    return "root_pass"


def continue_low(E,first,left,right,G,Q,k):
    """Use the previous actual cofactor coefficient to form the new irrational part."""
    out=[]
    for X,Y in first:
        internal=E["cross"](X,Y,k)
        for _,ax in E["shifts"](scale(G,left)):
            for budget in (Z,I):
                ay=sub(sub(budget,internal),ax)
                if ay[1]!=Q*right[1] or E["sign"](ay)<0 or E["sign"](sub(I,ay))<0:continue
                out.append((X+(ax,),Y+(ay,)))
    return out


def recover(E,p,G,A,k):
    x,y=sub(A,E["base"](p,G,k))
    assert not y
    return G+(x,)


def reverse_budgets(E,p,known):
    # Run the reverse rows with binary convolution values as the primitive branching.
    for X,Y in known:
        states=[((1,),(1,),(I,),(I,))]
        for k in range(1,10):
            states=E["stage"](p,states,k,force_coeff=Z if k==3 else None,against=(X[k],Y[k]))
            if not states:break
        if states:return True
    return False


def solve(job):
    (L,M,T,D),count3=job
    E=env(T,D)
    P=(I,z(L),z(M),A)
    Pb=tuple(map(E["conj"],P))
    S=ordinary_moments(P,E)
    assert S is not None
    states=[((1,),(1,),(I,),(I,))]
    for k in range(1,7):
        states=E["stage"](P,states,k,force_g=0 if k==3 else None,
                          force_coeff=Z if k==6 else None)
    counts=Counter(row3_states=count3)
    if not states:return counts,[]
    assert len(states)==1
    counts["row6_cores"]+=1
    G0,Q0,X,Y=states[0]
    G0,Q0=G0[:4],Q0[:4]
    result=[]
    U4=E["moment"](P,S,4)
    for c in E["bounds"](U4,4):
        P4,P4b=P+(z(c),),Pb+(z(c),)
        G4,Q4=recover(E,P4,G0,X[4],4),recover(E,P4b,Q0,Y[4],4)
        live7=continue_low(E,[(X,Y)],A,E["conj"](A),G4[-1],Q4[-1],7)
        if not live7:continue
        counts["row7_states"]+=len(live7)
        S4=S+[add(U4,z(-4*c))]
        U5=E["moment"](P4,S4,5)
        for d in E["bounds"](U5,5):
            P5,P5b=P4+(z(d),),P4b+(z(d),)
            G5,Q5=recover(E,P5,G4,X[5],5),recover(E,P5b,Q4,Y[5],5)
            live8=continue_low(E,live7,A,E["conj"](A),G5[-1],Q5[-1],8)
            if not live8:continue
            counts["row8_states"]+=len(live8)
            S5=S4+[add(U5,z(-5*d))]
            U6=E["moment"](P5,S5,6)
            for e in E["bounds"](U6,6):
                P6,P6b=P5+(z(e),),P5b+(z(e),)
                G6,Q6=recover(E,P6,G5,X[6],6),recover(E,P6b,Q5,Y[6],6)
                live9=continue_low(E,live8,A,E["conj"](A),G6[-1],Q6[-1],9)
                if not live9:continue
                counts["row9_states"]+=len(live9)
                S6=S5+[add(U6,z(-6*e))]
                U7=E["moment"](P6,S6,7)
                for g in E["bounds"](U7,7):
                    P7=P6+(z(g),)
                    S7=S6+[add(U7,z(-7*g))]
                    U8=E["moment"](P7,S7,8)
                    for h in E["bounds"](U8,8):
                        if abs(h*h-2*g)>17 or abs(-h**3+3*h*g-3*e)>33:continue
                        P9=P7+(z(h),I)
                        S8=S7+[add(U8,z(-8*h))]
                        if not E["house"](add(E["moment"](P9,S8,9),z(-9)),9):continue
                        rev=P9[::-1]
                        if ordinary_moments(rev,E) is None:continue
                        counts["newton_both"]+=1
                        if not reverse_budgets(E,rev,live9):continue
                        counts["coupled_mothers"]+=1
                        code=root_code(P9,E)
                        counts[code]+=1
                        result.append(((L,M,T,D,c,d,e,g,h),code))
    return counts,result


def independent_tree(record):
    row,tree,status=record
    assert status=="forward_dead"
    L,M,T,D,c,d,e,g,h=row
    E=env(T,D)
    P=(I,z(L),z(M),A,z(c),z(d),z(e),z(g),z(h),I)
    states=[((1,),(1,),(I,),(I,))]
    sizes=[]
    for k in range(1,tree["forward"]["first_empty"]+1):
        states=E["stage"](P,states,k,force_g=0 if k==3 else None,force_coeff=Z if k==6 else None)
        sizes.append(len(states))
    assert sizes==tree["forward"]["layer_sizes"]
    assert sizes[-1]==0
    return len(sizes)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("entry",type=Path)
    ap.add_argument("trees",type=Path)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=8)
    args=ap.parse_args()
    assert not args.out.exists()
    roots=list(core_rows())
    counts=Counter()
    result=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for ct,output in pool.map(solve,roots,chunksize=1):
            counts.update(ct)
            result.extend(output)
        data=b"".join(json.dumps(x,separators=(",",":")).encode()+b"\n" for x in sorted(result))
        reference=args.entry.read_bytes()
        assert data==reference, (len(result),counts)
        tree_raw=args.trees.read_bytes()
        survivors=[row for row,code in map(json.loads,reference.splitlines()) if code=="root_pass"]
        tree_records=list(map(json.loads,tree_raw.splitlines()))
        assert [r[0] for r in tree_records]==survivors
        deaths=Counter(pool.map(independent_tree,tree_records,chunksize=3))
    primary=json.loads(args.entry.with_name("r3-full.json").read_bytes())
    assert len(roots)==primary["input_cores"]
    assert all(counts[k]==v for k,v in primary["counts"].items())
    assert counts["row3_states"]==721
    output=dict(status="verified",first_cores=len(roots),counts=counts,independent_tree_deaths=deaths,
                mother_sha256=sha256(data).hexdigest(),tree_sha256=sha256(tree_raw).hexdigest())
    args.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    print(json.dumps(output,sort_keys=True))


if __name__=="__main__":
    main()
