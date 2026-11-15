"""Exact normalization audit of all stable promotions from 3x3 one-peak labels."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
import json
from exact import grid, peaks


def promote(labels,v,t):
    k=labels[v]
    assert k<t
    out=[x-1 if k<x<=t else x for x in labels]
    out[v]=t
    return tuple(out)


def bound(N,delta,F):
    return sum((Fraction(max(0,F-delta*q),N-q-1) for q in range(1,N-1)),Fraction(0))


def main():
    n=3;N=n*n;xy,adj=grid(n,n)
    counts=Counter();loads=Counter();source_min=Counter();triple_dist=Counter()
    for L in permutations(range(1,N+1)):
        pp=peaks(L,adj);counts[len(pp)]+=1
        if len(pp)!=1:continue
        u=pp[0]
        for v in range(N):
            if v==u or adj[u]>>v&1:continue
            M=max(L[x] for x in range(N) if adj[v]>>x&1)
            for t in range(M,N):
                image=promote(L,v,t);p=peaks(image,adj)
                assert set(p)=={u,v} and image[u]==N and image[v]==t
                loads[image]+=1
                triple_dist[abs(xy[u][0]-xy[v][0])+abs(xy[u][1]-xy[v][1])]+=1
    mass=Fraction(0);far_mass={1:Fraction(0),2:Fraction(0),3:Fraction(0)}
    max_weight=Fraction(0);hist=Counter()
    for image,c in loads.items():
        pp=peaks(image,adj);v=min(pp,key=lambda x:image[x]);t=image[v]
        assert c<=t-1
        w=Fraction(c,t-1);assert w<=1;mass+=w;max_weight=max(max_weight,w)
        u=next(x for x in pp if x!=v);d=abs(xy[u][0]-xy[v][0])+abs(xy[u][1]-xy[v][1]);hist[d]+=1
        for r in far_mass:
            if d>r:far_mass[r]+=w
    inverse_checks=0
    for image in permutations(range(1,N+1)):
        pp=peaks(image,adj)
        if len(pp)!=2:continue
        u=max(pp,key=lambda x:image[x]);v=min(pp,key=lambda x:image[x]);t=image[v]
        H=((1<<N)-1)^(1<<v)
        exposed=[x for x in peaks(image,adj,H) if x!=u]
        assert all(adj[v]>>x&1 for x in exposed)
        M=max(image[x] for x in range(N) if adj[v]>>x&1)
        B=max((image[x] for x in exposed),default=0)
        assert M-B==loads.get(image,0)
        inverse_checks+=1
    assert inverse_checks==counts[2]
    output={'inverse_formula_checks':inverse_checks,'n':n,'A1':counts[1],'A2':counts[2],'relation_triples':sum(loads.values()),
            'max_unweighted_inverse':max(loads.values()),'max_normalized_output_weight':str(max_weight),
            'total_weight':str(mass),'total_weight_divided_by_A1':str(mass/counts[1]),
            'far_weight':{str(k):str(v) for k,v in far_mass.items()},'distinct_output_distance_histogram':dict(hist),
            'weighted_ratio_bound':str(bound(N,4,N-1))}
    # Certified logarithm coefficient: f=7/8, Delta=4; log(1-x) via rational series.
    f=Fraction(7,8);x=f/4;K=50
    log_upper=-sum((x**k/k for k in range(1,K+1)),Fraction(0))
    log_lower=log_upper-x**(K+1)/((K+1)*(1-x))
    lower=f+(4-f)*log_lower;upper=f+(4-f)*log_upper
    assert Fraction(1,10)<lower<upper
    output['coefficient_epsilon_one_quarter']={'lower':str(lower),'upper':str(upper),
                                             'decimal_for_display':float((lower+upper)/2),
                                             'certified_greater_than':'1/10','series_terms':K}
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
