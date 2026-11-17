#!/usr/bin/env python3
"""Exact finite detailed-balance audit of the rank-biased prefix proposal.

The target is uniform on ordered valid prefixes with births exactly at ranks
1 and 2.  A proposed invalid prefix is a rejection/self-loop, precisely as
in the Metal kernel.  All arithmetic below is Fraction arithmetic.
"""
from __future__ import annotations

import itertools, json
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent

def neighbours(n, v):
    x,y=v%n,v//n
    if x: yield v-1
    if x+1<n: yield v+1
    if y: yield v-n
    if y+1<n: yield v+n

def valid_tau2(n, word):
    seen=set(); births=0
    for q,v in enumerate(word):
        if v in seen: return False
        born=not any(w in seen for w in neighbours(n,v))
        want=(q==0 or q==1)
        if born != want: return False
        births+=born; seen.add(v)
    return births==(1 if len(word)<2 else 2)

def states(n,r):
    return tuple(w for w in itertools.permutations(range(n*n),r) if valid_tau2(n,w))

def rankbiased(r):
    if r==1: return (Fraction(1),)
    return tuple(Fraction(1,4)+Fraction(1,2*r) if i<2 else Fraction(1,2*r) for i in range(r))

def proposed(word,i,v):
    y=list(word)
    try: j=y.index(v)
    except ValueError: j=len(y)
    oldi=y[i]
    y[i]=v
    if j<len(y): y[j]=oldi
    return tuple(y)

def kernel(n,r,p):
    ss=states(n,r); target=set(ss); N=n*n
    q={x:{} for x in ss}; accepted=set()
    for x in ss:
        for i,pi in enumerate(p):
            for v in range(N):
                y=proposed(x,i,v)
                if y in target:
                    accepted.add(y); q[x][y]=q[x].get(y,Fraction())+pi/N
                else:
                    q[x][x]=q[x].get(x,Fraction())+pi/N
        assert sum(q[x].values(),Fraction())==1
        assert set(q[x])<=target
    # Checking every nonzero entry suffices: if q(y,x)>0 but q(x,y)=0,
    # that same discrepancy is encountered when the loop reaches y.
    for x in ss:
        for y,mass in q[x].items():
            assert mass==q[y].get(x,Fraction()), (n,r,x,y,mass,q[y].get(x))
    # The accepted candidate set is a subset of the fixed target by design;
    # its support is not claimed irreducible or mixed.
    return {'states':len(ss),'nonzero_transition_pairs':sum(len(z) for z in q.values()),
            'accepted_target_states':len(accepted),'all_rows_sum_one':True,
            'offdiagonal_detailed_balance':True,'proposal_never_leaves_target':True}

def main():
    out={'target':'uniform ordered tau=2 valid prefix','exact_arithmetic':'fractions.Fraction',
         'rankbiased_rule':'r=1: p0=1; r>=2: p0=p1=1/4+1/(2r), pi=1/(2r) for i>=2',
         'checks':[]}
    # n=2 is complete at every meaningful prefix rank; n=3 is deliberately
    # restricted to small r, where full enumerations remain independent and
    # transparent rather than a mixing experiment.
    for n,rs in ((2,range(1,5)),(3,range(1,6))):
        for r in rs:
            p=rankbiased(r)
            a=kernel(n,r,p)
            # A second arbitrary fixed weight vector checks the general
            # fixed-rank-weight symmetry, beyond the particular r025 choice.
            generic=tuple(Fraction(i+1,sum(range(1,r+1))) for i in range(r))
            b=kernel(n,r,generic)
            out['checks'].append({'n':n,'r':r,'rankbiased_weights':[str(x) for x in p],
                                  'rankbiased':a,'generic_weights':[str(x) for x in generic],
                                  'generic_fixed_weight':b})
    (HERE/'exact-prefix-balance.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
