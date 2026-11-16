#!/usr/bin/env python3
"""Independent Fraction implementation of the n=3 first-merge decomposition."""
import json
from fractions import Fraction
from math import comb
from pathlib import Path

n=3; N=n*n; cap=1<<N; full=cap-1; L=1<<18
adj=[0]*N
for y in range(n):
    for x in range(n):
        v=y*n+x
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            X,Y=x+dx,y+dy
            if 0<=X<n and 0<=Y<n: adj[v]|=1<<(Y*n+X)

def comp(seed, allowed):
    seen=front=seed
    while front:
        b=front&-front;front-=b
        extra=adj[b.bit_length()-1]&allowed&~seen
        seen|=extra;front|=extra
    return seen

a=[0]*cap; mx=[0]*cap; my=[0]*cap; m2=[0]*cap
for v in range(N):
    x,y=v%n,v//n; a[1<<v]=1;mx[1<<v]=x;my[1<<v]=y;m2[1<<v]=x*x+y*y
for S in range(1,cap):
    if S&(S-1)==0: continue
    for v in range(N):
        if S>>v&1:
            T=S^(1<<v)
            if adj[v]&T:
                a[S]+=a[T];mx[S]+=mx[T];my[S]+=my[T];m2[S]+=m2[T]
h=[0]*cap;h[full]=1
for S in range(full-1,0,-1):
    for v in range(N):
        if not(S>>v&1) and adj[v]&S:h[S]+=h[S|(1<<v)]

W=D=R2=0; mean=Fraction(); within=Fraction(); between=Fraction()
for S in range(1,cap):
    if not a[S] or not h[S]:continue
    s=S.bit_count()
    for w in range(N):
        if not(S>>w&1):continue
        T=S^(1<<w); A=comp(T&-T,T);B=T^A
        if not B or comp(B&-B,B)!=B:continue
        aa,bb=a[A],a[B]; base=h[S]*comb(s-1,A.bit_count()); weight=base*aa*bb
        rr=base*(bb*m2[A]+aa*m2[B])
        dd=base*(bb*m2[A]+aa*m2[B]-2*(mx[A]*mx[B]+my[A]*my[B]))
        ma=Fraction(mx[A]*mx[A]+my[A]*my[A],aa*aa)
        mb=Fraction(mx[B]*mx[B]+my[B]*my[B],bb*bb)
        va=Fraction(m2[A],aa)-ma; vb=Fraction(m2[B],bb)-mb
        W+=weight;D+=dd;R2+=rr;mean+=weight*(ma+mb);within+=weight*(va+vb);between+=dd-weight*(va+vb)

c=json.loads(Path('mean-variance-n3.json').read_text())['totals']
def interval(lo,hi):
    return Fraction(int(lo),L*L),Fraction(int(hi),L*L)
meanlo,meanhi=interval(c['mean_square_lower_numerator'],c['mean_square_upper_numerator'])
winlo,winhi=interval(c['within_lower_numerator'],c['within_upper_numerator'])
betlo,bethi=interval(c['between_lower_numerator'],c['between_upper_numerator'])
assert (W,D,R2)==(int(c['W']),int(c['D']),int(c['root_second_moment']))
assert meanlo<=mean<=meanhi and winlo<=within<=winhi and betlo<=between<=bethi
assert within+between==D
out={'status':'passed','n':3,'W':W,'D':D,'exact_mean_square':str(mean),'exact_within':str(within),'exact_between':str(between),'fixed_point_intervals_contain_exact':True}
Path('audit-mean-variance-n3.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
