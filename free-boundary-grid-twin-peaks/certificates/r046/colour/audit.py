"""Exact small gates and uniform finite-n4 sampling for the colour identity."""
from collections import Counter
from itertools import combinations, permutations, product
import json
import math
from pathlib import Path
import random
import sys

def geometry(n):
    vertices=range(n*n)
    nbr=[tuple(w for w in vertices if abs(v%n-w%n)+abs(v//n-w//n)==1) for v in vertices]
    masks=[sum(1<<w for w in ns) for ns in nbr]
    classes=[[v for v in vertices if (v%n+2*(v//n))%5==c] for c in range(5)]
    decomposition=[]
    for centres in classes:
        patches=[(v,*nbr[v]) for v in centres]
        occupied={p for patch in patches for p in patch}
        assert len(occupied)==sum(len(patch) for patch in patches)
        decomposition.append((centres,tuple(v for v in vertices if v not in occupied)))
    return nbr,masks,decomposition

def peaks(x,nbr):
    return tuple(v for v,ns in enumerate(nbr) if all(x[v]>x[w] for w in ns))

def coefficients(x,v,nbr,one):
    maximum=max(x[w] for w in nbr[v])
    exposed=sorted((x[w] for w in nbr[v] if all(x[w]>x[z] for z in nbr[w] if z!=v)),reverse=True)
    first=exposed[0] if exposed else 0
    second=exposed[1] if len(exposed)>1 else 0
    weights=[maximum-first,one-maximum+first-second]
    weights.extend(b-(exposed[j+1] if j+1<len(exposed) else 0)
                   for j,b in enumerate(exposed) if j>=1)
    return weights

def low_coefficients(x,nbr,centres,fixed,one):
    local=[coefficients(x,v,nbr,one) for v in centres]
    fixed_peaks=sum(all(x[v]>x[u] for u in nbr[v]) for v in fixed)
    minimums=[int(not row[0]) for row in local]
    m=fixed_peaks+sum(minimums)
    extras=[1,0,0]
    for j,row in zip(minimums,local):
        shifted=row[j:]+[0,0]
        a,b=shifted[1]/shifted[0],shifted[2]/shifted[0]
        extras=[1,extras[1]+a,extras[2]+b+extras[1]*a]
    return m,extras

def direct_gate(rng):
    n=3
    nbr,_,parts=geometry(n)
    checked=0
    for _ in range(140):
        ranks=list(range(1,n*n+1))
        rng.shuffle(ranks)
        for centres,fixed in parts:
            z=sorted({0,n*n+1}|{ranks[v] for v in range(n*n) if v not in centres})
            intervals=[((a+b)/2,b-a) for a,b in zip(z,z[1:])]
            actual=Counter()
            for config in product(intervals,repeat=len(centres)):
                y=ranks[:]
                weight=1
                for v,(x,length) in zip(centres,config):
                    y[v]=x
                    weight*=length
                actual[len(peaks(y,nbr))]+=weight
            base=sum(all(ranks[v]>ranks[u] for u in nbr[v]) for v in fixed)
            predicted=Counter({base:1})
            for v in centres:
                nxt=Counter()
                for a,w in predicted.items():
                    for j,p in enumerate(coefficients(ranks,v,nbr,n*n+1)):
                        nxt[a+j]+=w*p
                predicted=nxt
            assert actual==predicted
            checked+=1
    return checked

def exhaustive_slack():
    n=3
    nbr,_,parts=geometry(n)
    slack=Counter()
    by_distance=Counter()
    for ranks in permutations(range(1,10)):
        roots=peaks(ranks,nbr)
        if len(roots)!=2: continue
        this=[]
        for centres,fixed in parts:
            m=sum(all(ranks[v]>ranks[u] for u in nbr[v]) for v in fixed)
            m+=sum(coefficients(ranks,v,nbr,10)[0]==0 for v in centres)
            this.append(m)
        assert max(this)<=2
        u,v=roots
        d=abs(u%n-v%n)+abs(u//n-v//n)
        by_distance[d]+=1
        slack[d,min(this)]+=1
    assert sum(by_distance.values())==143808
    return {str(d):{str(m):slack[d,m] for m in range(3)} for d in sorted(by_distance)}

def sample_n4(rng,draws):
    n=4
    N=n*n
    nbr,nmask,parts=geometry(n)
    full=(1<<N)-1
    completion=[(0,0,0)]*(full+1)
    completion[full]=(1,0,0)
    for s in range(full-1,-1,-1):
        row=[0,0,0]
        remainder=full^s
        while remainder:
            b=remainder&-remainder
            remainder^=b
            p=b.bit_length()-1
            child=completion[s|b]
            birth=int(not(nmask[p]&s))
            for k in range(birth,3): row[k]+=child[k-birth]
        completion[s]=tuple(row)
    assert completion[0][2]==779299174080
    minima=Counter()
    selected=[]
    for _ in range(draws):
        values=sorted((rng.random() for _ in range(N)),reverse=True)
        labels=[0]*N
        s,remaining_peaks=0,2
        for value in values:
            choice=rng.randrange(completion[s][remaining_peaks])
            remainder=full^s
            while remainder:
                b=remainder&-remainder
                remainder^=b
                p=b.bit_length()-1
                birth=int(not(nmask[p]&s))
                z=0 if birth>remaining_peaks else completion[s|b][remaining_peaks-birth]
                if choice<z: break
                choice-=z
            labels[p]=value
            s|=b
            remaining_peaks-=birth
        assert s==full and remaining_peaks==0 and len(peaks(labels,nbr))==2
        stats=[low_coefficients(labels,nbr,c,f,1) for c,f in parts]
        minima[min(m for m,w in stats)]+=1
        active=[w[2-m]/N**(2-m) for m,w in stats if m<2]
        if active:selected.append(min(active))
    selected.sort()
    return dict(draws=draws,normalizer=completion[0][2],minimum_m_counts=minima,
                scaled_weight_quantiles={str(q):selected[int(q*(len(selected)-1))] for q in (.25,.5,.75,.95,.99)})

if __name__=='__main__':
    rng=random.Random(20261118046)
    output=dict(direct_skeletons=direct_gate(rng),n3_exact_min_m_by_distance=exhaustive_slack(),
                n4_numerical=sample_n4(rng,16384),seed=20261118046,python=sys.version)
    Path(sys.argv[1]).write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output))
