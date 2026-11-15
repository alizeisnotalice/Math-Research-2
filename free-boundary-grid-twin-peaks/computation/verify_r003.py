"""Independent checks of round 3 certificates and the all-n gap-one family."""
from collections import deque
from fractions import Fraction
from pathlib import Path
import json
from exact import grid, peaks

BASE=Path(__file__).parent/'results'
def gap_one(n):
    N=n*n;_,adj=grid(n,n);full=(1<<N)-1
    v=N-1;x=N-2;y=N-1-n;q=N-2-n;r=N-1-2*n
    A=full^sum(1<<w for w in (v,x,y,q,r));todo=deque([0]);seen={0};order=[]
    while todo:
        w=todo.popleft();order.append(w)
        for z in range(N):
            if (adj[w]>>z&1) and (A>>z&1) and z not in seen:seen.add(z);todo.append(z)
    assert len(order)==N-5 and not(adj[v]&A)
    total=Fraction(0)
    for j in range(1,len(order)+1):
        seq=order[:j]+[v]+order[j:]+[x,y,q,r];L=[0]*N
        for k,w in enumerate(seq):L[w]=N-k
        assert peaks(L,adj)==[0,v]
        exposed=[z for z in peaks(L,adj,full^(1<<v)) if z!=0]
        M=max(L[z] for z in range(N) if adj[v]>>z&1);b=max(L[z] for z in exposed)
        assert M==L[x] and b==L[y] and M-b==1
        total+=Fraction(M-b,L[v]-1)
    assert total==sum((Fraction(1,k) for k in range(4,N-1)),Fraction(0))
    return {'n':n,'checked_insertions':N-5,'distance':2*(n-1),'beta_mean':str(total/(N-5))}

A2={3:143808,4:779299174080,5:7063261481362958280000}
W={3:Fraction(262272,5),4:Fraction(5570097564160,39),5:Fraction(28875173046780972423711380,46189)}
expected={3:79360,4:266753826624,5:1353615021505617063360}
rows=[]
for n in (3,4,5):
    z=json.loads((BASE/f'singleton-n{n}-r003.json').read_text())
    assert z['beta_positive_count']==expected[n] and 2*z['both_singleton_count']==z['A1']
    assert 0<W[n]<=expected[n]<A2[n]
    rows.append({'n':n,'support_probability':str(Fraction(expected[n],A2[n])),
                 'mean_beta_given_positive':str(W[n]/expected[n])})
z=json.loads((BASE/'shortest-lift-r003.json').read_text())
assert z['checked_A2']==A2[3] and sum(z['distance_counts'].values())==A2[3]
for d,c in z['max_fiber_by_distance'].items():assert c<=z['universal_fiber_bounds'][d]
walk=json.loads((BASE/'beta-walks-r003.json').read_text())
for z in walk['results']:
    n=z['n'];L=z['best_labels'];assert sorted(L)==list(range(1,n*n+1))
    _,adj=grid(n,n);p=peaks(L,adj);assert len(p)==2
    d=abs(p[0]%n-p[1]%n)+abs(p[0]//n-p[1]//n)
    assert d==z['best_distance'] and sum(z['distance_histogram'])==z['samples']
    assert all(a+b==c for a,b,c in zip(z['first_half_histogram'],z['second_half_histogram'],z['distance_histogram']))
    assert 0<=z['beta_positive_samples']<=z['samples']
assert sum(z['steps'] for z in walk['results'])==280000000
print(json.dumps({'verified':True,'support_exact_ratios':rows,'gap_one_fibers':[gap_one(n) for n in (5,7,11)],
                  'walks_checked':len(walk['results']),'proposals':280000000},indent=2))
