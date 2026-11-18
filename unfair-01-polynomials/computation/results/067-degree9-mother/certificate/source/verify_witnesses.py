#!/usr/bin/env python3
"""Check the two short prefix witnesses in independent integral coordinates."""

import argparse
from collections import Counter
import json
from pathlib import Path

from singleton_probe import ONE,ZERO,ci,conj,extend,newton_prefix,plus,root_code,scale


def direct(P,G):
    output=[]
    for k in range(len(G)):
        u,v=0,0
        for j in range(min(k,len(P)-1)+1):
            w,z=P[j]
            u,v=u+G[k-j]*w,v+G[k-j]*z
        output.append((u,v))
    return output


def run_tree(P,r,forward,depth=80,cap=10000):
    rows=[((1,),(1,),(ONE,),(ONE,))]
    trace=[]
    Db=5
    Pb=list(map(conj,P))
    for k in range(1,depth+1):
        next_rows=[]
        for G,Q,A,B in rows:
            x=ZERO
            y=ZERO
            for j in range(1,min(9,k)+1):
                x=plus(x,scale(G[k-j],P[j]))
                y=plus(y,scale(Q[k-j],Pb[j]))
            for u,v,a,b in extend(A,B,x,y,k,Db):
                if k==(9-r if forward else r) and (a!=ZERO or b!=ZERO):continue
                if forward and k==9-2*r and (u or v):continue
                next_rows.append((G+(u,),Q+(v,),A+(a,),B+(b,)))
        trace.append(len(next_rows))
        if not next_rows:return dict(status="dead",first_empty=k,layer_sizes=trace)
        if len(next_rows)>cap:return dict(status="unknown_cap",last_layer=k,layer_sizes=trace)
        rows=next_rows
    return dict(status="unknown_depth",last_layer=depth,layer_sizes=trace)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("source",type=Path)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    assert not args.out.exists()
    source=json.loads(args.source.read_bytes())
    report={}
    for r in (2,4):
        w=source[str(r)]["witness"]
        assert w["D"]==5
        P=list(map(tuple,w["P"]))
        assert P==[ONE]+[(1,1) if k==r else ZERO for k in range(1,9)]+[ONE]
        G,Q,A,B,U,V,C,F=w["simultaneous_prefix"]
        A,B,C,F=[list(map(tuple,a)) for a in (A,B,C,F)]
        assert all(len(a)==10 for a in (G,Q,A,B,U,V,C,F))
        assert direct(P,G)==A and direct(list(map(conj,P)),Q)==B
        assert direct(P[::-1],U)==C and direct(list(map(conj,P[::-1])),V)==F
        # For D=5 use the full integral basis beta=(1+sqrt(5))/2.
        def integral(x):
            u,v=x
            assert (u-v)%2==0
            return (u-v)//2,v
        def mul(x,y):
            m,n=x
            p,q=y
            return m*p+n*q,m*q+n*p+n*q
        def nonnegative(z):
            m,n=z
            u=2*m+n
            if u>=0 and n>=0:return True
            if u<=0 and n<=0:return not(u or n)
            return u*u>5*n*n if u>0 else u*u<5*n*n
        def products(X,Y):
            output=[]
            x,y=list(map(integral,X)),list(map(integral,Y))
            for k in range(10):
                terms=[mul(x[j],y[k-j]) for j in range(k+1)]
                output.append(tuple(sum(t[i] for t in terms) for i in (0,1)))
            return output
        for row in (A,B,C,F):
            assert all(nonnegative(z) and nonnegative((1-z[0],-z[1])) for z in map(integral,row))
        assert all(c in ((0,0),(1,0)) for c in products(A,B)+products(C,F))
        assert all(not ((A[k]!=ZERO and F[k]!=ZERO) or (B[k]!=ZERO and C[k]!=ZERO)) for k in range(1,10))
        assert G[9-2*r]==Q[9-2*r]==0
        assert A[9-r]==B[9-r]==C[r]==F[r]==ZERO
        assert newton_prefix(P,5) and newton_prefix(P[::-1],5)
        assert root_code(P,5)=="root_pass"
        # For all positive x, the only negative coefficient is beta' in (-1,0).
        assert (2,0)==P[0]==P[9]
        report[r]=dict(scope="both depth-nine prefixes, exact root annulus and full positive axis",
                       forward_tree=run_tree(P,r,True),reverse_tree=run_tree(P[::-1],r,False))
    args.out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print(json.dumps(report,sort_keys=True))


if __name__=="__main__":
    main()
