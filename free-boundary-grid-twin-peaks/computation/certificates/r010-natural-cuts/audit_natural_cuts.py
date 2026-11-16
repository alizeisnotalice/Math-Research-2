#!/usr/bin/env python3
"""Exact n<=3 audit of the natural-forest cut lower bound.

For a left/right cut A|B, roots are permitted only at *safe* outer vertices:
their sole inward neighbour is in the same side.  This qualification is
essential; merely placing a root in the outermost column does not prevent a
cross-edge from destroying that global peak.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import permutations, product
import json
from math import comb, factorial
from pathlib import Path


def graph(n):
    V = range(n*n)
    adj = [0]*(n*n)
    for y in range(n):
        for x in range(n):
            v=y*n+x
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                X,Y=x+dx,y+dy
                if 0<=X<n and 0<=Y<n: adj[v] |= 1 << (Y*n+X)
    return adj


def connected(mask,adj):
    if not mask: return False
    seen=mask & -mask
    frontier=seen
    while frontier:
        nb=0
        x=frontier
        while x:
            b=x & -x; x-=b; nb |= adj[b.bit_length()-1]
        frontier=nb & mask & ~seen
        seen |= frontier
    return seen==mask


def onepeak_count(mask,allowed,adj):
    """Exact count of label orders with one peak in allowed."""
    memo={}
    def f(t):
        if t in memo: return memo[t]
        if t & (t-1)==0:
            z=int(bool(t & allowed)); memo[t]=z; return z
        z=0; x=t
        while x:
            b=x & -x; x-=b
            if connected(t^b,adj): z += f(t^b)
        memo[t]=z
        return z
    return f(mask)


def safe_left(mask,n):
    return sum((1<<(y*n)) for y in range(n)
               if mask&(1<<(y*n)) and mask&(1<<(y*n+1)))


def safe_right(mask,n):
    return sum((1<<(y*n+n-1)) for y in range(n)
               if mask&(1<<(y*n+n-1)) and mask&(1<<(y*n+n-2)))


def safe_left_all(mask,n,adj):
    return sum((1<<(y*n)) for y in range(n)
               if mask&(1<<(y*n)) and not (adj[y*n] & ~mask))


def safe_right_all(mask,n,adj):
    return sum((1<<(y*n+n-1)) for y in range(n)
               if mask&(1<<(y*n+n-1)) and not (adj[y*n+n-1] & ~mask))


def cross_edges(A,adj):
    B=((1<<len(adj))-1)^A
    return sum((adj[v]&B).bit_count() for v in range(len(adj)) if A>>v&1)


def peakset(lab,adj):
    return tuple(v for v in range(len(adj))
                 if all(lab[v]>lab[u] for u in range(len(adj)) if adj[v]>>u&1))


def induced_peakset(lab,mask,adj):
    return tuple(v for v in range(len(adj)) if mask>>v&1
                 if all(lab[v]>lab[u] for u in range(len(adj))
                        if (mask>>u&1) and (adj[v]>>u&1)))


def forest_cut_prob(lab,A,adj):
    """Conditional probability that every higher-neighbour arrow stays in side."""
    roots=peakset(lab,adj)
    if len(roots)!=2: return Fraction(0)
    B=((1<<len(adj))-1)^A
    p=Fraction(1)
    for v in range(len(adj)):
        higher=[u for u in range(len(adj)) if adj[v]>>u&1 and lab[u]>lab[v]]
        if not higher: continue
        inside=[u for u in higher if ((A>>u)&1)==((A>>v)&1)]
        if not inside: return Fraction(0)
        p*=Fraction(len(inside),len(higher))
    return p


def audit(n):
    N=n*n; full=(1<<N)-1; adj=graph(n)
    left=sum(1<<(y*n) for y in range(n)); right=sum(1<<(y*n+n-1) for y in range(n))
    cuts=[]
    for A in range(1<<N):
        if A&left!=left or A&right: continue
        B=full^A
        if not connected(A,adj) or not connected(B,adj): continue
        FL=safe_left(A,n); FR=safe_right(B,n)
        if not FL or not FR: continue
        a=onepeak_count(A,FL,adj); b=onepeak_count(B,FR,adj)
        if a and b:
            c=comb(N,A.bit_count())*a*b
            cuts.append((A,c,cross_edges(A,adj),FL,FR))
    weighted=sum(Fraction(c,2**e) for _,c,e,_,_ in cuts)
    # Full-label audit only n<=3.
    event_counts=defaultdict(int); forest_mass=Fraction(); far=0; k2=0
    for lab in permutations(range(N)):
        P=peakset(lab,adj)
        if len(P)!=2: continue
        k2+=1
        rho=abs(P[0]%n-P[1]%n)+abs(P[0]//n-P[1]//n)
        if rho>=n-1: far+=1
        for A,c,e,FL,FR in cuts:
            B=full^A
            # The constructed event: each side has exactly one internal peak in its safe set.
            PA=induced_peakset(lab,A,adj); PB=induced_peakset(lab,B,adj)
            if len(PA)==len(PB)==1 and (FL>>PA[0]&1) and (FR>>PB[0]&1):
                # Safe roots have no cross neighbour, hence these are exactly
                # the two global peaks of this full labelling.
                assert set(P)==set(PA+PB)
                event_counts[A]+=1
                theta=forest_cut_prob(lab,A,adj)
                assert theta >= Fraction(1,2**e), (n,A,lab,theta,e)
                forest_mass += theta
    assert all(event_counts[A]==c for A,c,_,_,_ in cuts), (n,cuts,event_counts)
    assert forest_mass <= far <= k2
    assert weighted <= forest_mass
    return {"n":n,"A2":k2,"far_rho_ge_n_minus_1":far,"cuts":[
        {"A_mask":A,"raw_count":c,"cross_edges":e,"safe_left":FL,"safe_right":FR}
        for A,c,e,FL,FR in cuts],"weighted_lower":[weighted.numerator,weighted.denominator],
        "forest_mass":[forest_mass.numerator,forest_mass.denominator],"status":"passed"}


def cuts_only(n):
    """Exact construction ledger, without enumerating full labellings."""
    N=n*n; full=(1<<N)-1; adj=graph(n)
    left=sum(1<<(y*n) for y in range(n)); right=sum(1<<(y*n+n-1) for y in range(n))
    rows=[]
    # Left and right columns are fixed; only n(n-2) bits are free.
    free=[v for v in range(N) if not (left>>v&1 or right>>v&1)]
    for code in range(1<<len(free)):
        A=left
        for i,v in enumerate(free):
            if code>>i&1:A|=1<<v
        B=full^A
        if not connected(A,adj) or not connected(B,adj):continue
        FL=safe_left(A,n);FR=safe_right(B,n)
        if not FL or not FR:continue
        a=onepeak_count(A,FL,adj);bb=onepeak_count(B,FR,adj)
        if not (a and bb):continue
        e=cross_edges(A,adj);c=comb(N,A.bit_count())*a*bb
        rows.append((A.bit_count(),e,c,A,FL,FR))
    by={}
    for size,e,c,*_ in rows:
        z=by.setdefault((size,e),0);by[(size,e)]=z+c
    weighted=sum(Fraction(c,2**e) for _,e,c,*_ in rows)
    return {"n":n,"cut_count":len(rows),"raw_sum":sum(x[2] for x in rows),
       "weighted_lower":[weighted.numerator,weighted.denominator],
       "by_size_cross":[{"A_size":k[0],"cross_edges":k[1],"raw_count":v,
                         "weighted_numerator_over_2e":v} for k,v in sorted(by.items())],
       "status":"exact_construction_ledger"}


def allsafe_audit(n):
    N=n*n; full=(1<<N)-1; adj=graph(n); cuts=[]
    for A in range(1,full):
        B=full^A
        if not connected(A,adj) or not connected(B,adj): continue
        FL=safe_left_all(A,n,adj); FR=safe_right_all(B,n,adj)
        if not FL or not FR: continue
        a=onepeak_count(A,FL,adj);b=onepeak_count(B,FR,adj)
        if a and b:cuts.append((A,comb(N,A.bit_count())*a*b,cross_edges(A,adj),FL,FR))
    weighted=sum(Fraction(c,2**e) for _,c,e,_,_ in cuts)
    counts=defaultdict(int);mass=Fraction();far=k2=0
    for lab in permutations(range(N)):
        P=peakset(lab,adj)
        if len(P)!=2:continue
        k2+=1;rho=abs(P[0]%n-P[1]%n)+abs(P[0]//n-P[1]//n)
        if rho>=n-1:far+=1
        per=Fraction()
        for A,c,e,FL,FR in cuts:
            B=full^A; PA=induced_peakset(lab,A,adj);PB=induced_peakset(lab,B,adj)
            if len(PA)==len(PB)==1 and (FL>>PA[0]&1) and (FR>>PB[0]&1):
                assert set(P)==set(PA+PB)
                counts[A]+=1;theta=forest_cut_prob(lab,A,adj);assert theta>=Fraction(1,2**e)
                per+=theta;mass+=theta
        assert per<=1
    assert all(counts[A]==c for A,c,_,_,_ in cuts)
    assert weighted<=mass<=far<=k2
    return {'n':n,'cuts':len(cuts),'raw_sum':sum(c for _,c,*_ in cuts),
      'weighted_lower':[weighted.numerator,weighted.denominator],
      'forest_mass':[mass.numerator,mass.denominator],'A2':k2,'far':far,'status':'passed'}


if __name__=='__main__':
    out={str(n):audit(n) for n in (2,3)}
    out.update({f"cuts_{n}":cuts_only(n) for n in (4,5)})
    out['allsafe_n3']=allsafe_audit(3)
    Path('audit-n2-n3.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
