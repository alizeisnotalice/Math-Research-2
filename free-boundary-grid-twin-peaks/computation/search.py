"""Biased constructive exploration. Its hit rates are NOT conditional probabilities."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import sys
import time
from exact import grid, vertices


def independent_peaks(labels,n):
    result=[]
    for y in range(n):
        for x in range(n):
            value=labels[y*n+x]
            neighbors=[labels[b*n+a] for a,b in ((x-1,y),(x+1,y),(x,y-1),(x,y+1))
                       if 0<=a<n and 0<=b<n]
            if all(value>other for other in neighbors):result.append(y*n+x)
    return result


def seeded_growth(n,seeds,rng,adj):
    # A construction law, not the uniform conditional law.
    order=list(seeds);used=set(seeds)
    while len(order)<n*n:
        boundary={u for v in used for u in vertices(adj[v]) if u not in used}
        v=rng.choice(sorted(boundary));used.add(v);order.append(v)
    labels=[0]*(n*n)
    for k,v in enumerate(order):labels[v]=n*n-k
    return labels


def describe(labels,n):
    pp=independent_peaks(labels,n)
    assert sorted(labels)==list(range(1,n*n+1))
    if len(pp)!=2:return {'peak_count':len(pp)}
    # Parent = largest labelled higher neighbor. Deterministic diagnostic forest.
    xy,adj=grid(n,n);roots={v:v for v in pp}
    for v in sorted(range(n*n),key=lambda v:-labels[v]):
        if v not in roots:
            parent=max((u for u in vertices(adj[v]) if labels[u]>labels[v]),key=lambda u:labels[u])
            roots[v]=roots[parent]
    boundary=sum(roots[u]!=roots[v] for v in range(n*n) for u in vertices(adj[v]) if u>v)
    return {'peak_count':2,'peaks':pp,
        'distance':sum(abs(a-b) for a,b in zip(xy[pp[0]],xy[pp[1]])),
        'basin_sizes':dict(Counter(roots.values())), 'deterministic_forest_cut_edges':boundary,
        'labels_row_major':labels}


def search(n,kind,seed,steps):
    rng=random.Random(seed);xy,bits=grid(n,n);adj=[list(vertices(a)) for a in bits];N=n*n
    if kind=='random':labels=rng.sample(range(1,N+1),N)
    else:
        seeds=(0,N-1) if kind=='far_growth' else (0,2)
        labels=seeded_growth(n,seeds,rng,bits)
    peak=lambda v:all(labels[v]>labels[u] for u in adj[v])
    pp={v for v in range(N) if peak(v)}
    def distance(p):
        if len(p)!=2:return 0
        u,v=p;return abs(u%n-v%n)+abs(u//n-v//n)
    def energy(p):return 10*abs(len(p)-2)-distance(p)/n
    best=labels[:] if len(pp)==2 else None;best_d=distance(pp) if best else -1
    E=energy(pp);accepted=hits=0;minK=len(pp);trace=[]
    start=time.perf_counter()
    for it in range(steps):
        v=rng.randrange(N)
        u=rng.choice(adj[v]) if rng.random()<.5 else rng.randrange(N)
        if u==v:continue
        affected={u,v,*adj[u],*adj[v]}
        labels[u],labels[v]=labels[v],labels[u]
        changed=(pp-affected)|{w for w in affected if peak(w)}
        en=energy(changed)
        beta=.2+7.8*(it/max(1,steps-1))
        if en<=E or rng.random()<math.exp(-beta*(en-E)):
            pp=changed;E=en;accepted+=1
        else:labels[u],labels[v]=labels[v],labels[u]
        minK=min(minK,len(pp))
        if len(pp)==2:
            hits+=1;d=distance(pp)
            if d>best_d:best=labels[:];best_d=d
        if (it+1)%(steps//10)==0:
            assert set(independent_peaks(labels,n))==pp
            trace.append({'step':it+1,'K':len(pp),'distance_if_two':distance(pp)})
    result={'n':n,'start_kind':kind,'seed':seed,'steps':steps,'accepted':accepted,
        'two_peak_visits_not_probability':hits,'minimum_peak_count':minK,
        'best':describe(best,n) if best else None,'terminal':describe(labels,n),
        'trace':trace,'wall_seconds':time.perf_counter()-start}
    if best:
        assert result['best']['peak_count']==2 and result['best']['distance']==best_d
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--steps',type=int,default=200000)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'purpose':'biased construction discovery only; no mixing or probability claims',
        'energy':'10*abs(K-2)-rho/n when K=2; 10*abs(K-2) otherwise',
        'schedule':'beta linearly increases from 0.2 to 8', 'runs':[]}
    for n in (8,12,16,24,32):
        for k,kind in enumerate(('random','near_growth','far_growth')):
            r=search(n,kind,2026111500+10*n+k,a.steps);out['runs'].append(r)
            a.output.write_text(json.dumps(out,indent=2)+'\n')
            print(json.dumps({key:r[key] for key in ('n','start_kind','minimum_peak_count','two_peak_visits_not_probability','wall_seconds')}),flush=True)
if __name__=='__main__':main()
