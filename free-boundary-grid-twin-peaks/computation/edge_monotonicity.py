#!/usr/bin/env python3
"""Exact small-graph test of internal-edge monotonicity of p."""
from fractions import Fraction
from functools import lru_cache
import json

def bits(x):
    while x:
        b=x&-x; x-=b; yield b

def connected(adj):
    n=len(adj); seen=q=1
    while q:
        b=q&-q; q-=b; v=b.bit_length()-1
        new=adj[v]&~seen; seen|=new; q|=new
    return seen==(1<<n)-1

def p(adj):
    n=len(adj); full=(1<<n)-1; root=1; F=adj[0]&~1
    @lru_cache(None)
    def h(T):
        if T==full:return 1
        available=0
        for v in range(n):
            if not (T>>v)&1 and adj[v]&T:available|=1<<v
        return sum(h(T|b) for b in bits(available))
    H=h(root); G=0
    for xb in bits(F):
        x=xb.bit_length()-1
        for yb in bits(adj[x]&~root&~F):G+=h(root|xb|yb)
    return Fraction(G,H)

def adj_of(n,code):
    a=[0]*n;k=0
    for i in range(n):
      for j in range(i+1,n):
        if(code>>k)&1:a[i]|=1<<j;a[j]|=1<<i
        k+=1
    return a

def main():
    rows=[]
    for n in range(4,7):
      cnt=0
      for code in range(1<<(n*(n-1)//2)):
        a=adj_of(n,code)
        if not connected(a) or not (1<=a[0].bit_count()<n-1):continue
        old=p(a)
        for i in range(1,n):
          for j in range(i+1,n):
            if (a[i]>>j)&1:continue
            b=a.copy();b[i]|=1<<j;b[j]|=1<<i
            new=p(b);cnt+=1
            if new<old:
              print(json.dumps({'status':'counterexample','n':n,'code':code,'edge':[i,j],
                 'adj':a,'p_before':[old.numerator,old.denominator],
                 'p_after':[new.numerator,new.denominator]},indent=2));return
      rows.append({'n':n,'comparisons':cnt})
    print(json.dumps({'status':'passed_through_n6','rows':rows},indent=2))
if __name__=='__main__':main()
