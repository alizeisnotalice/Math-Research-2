"""Independent finite result and D4 coverage checks for round 4."""
from pathlib import Path
from fractions import Fraction
from math import factorial
import json
B=Path(__file__).parent/'results'
A2={3:143808,5:7063261481362958280000,6:304351937651819767235871772704107520}
SUP={3:79360,5:1353615021505617063360,6:28664968745328849344136008673552640}
rows=[]
def transform(n,v,op):
    x,y=v%n,v//n
    if op>=4:x=n-1-x
    for _ in range(op%4):x,y=n-1-y,x
    return y*n+x
for n in (3,5,6):
    z=json.loads((B/f'beta-zero-orbits-n{n}-r004.json').read_text());edges={};total=0
    for o in z['orbits']:
        v,x=o['representative'];assert abs(v%n-x%n)+abs(v//n-x//n)==1
        orbit={(transform(n,v,k),transform(n,x,k)) for k in range(8)}
        assert len(orbit)==o['multiplicity'] and not(orbit&edges.keys())
        low=int(o['beta0']);assert low==int(o['E_K2'])-int(o['v_high'])>=0
        for e in orbit:edges[e]=low
        total+=len(orbit)*low
    assert len(edges)==4*n*(n-1)
    assert all(edges[v,x]==edges[x,v] for v,x in edges)
    assert total==int(z['beta0_total'])==A2[n]-SUP[n]
    rows.append({'n':n,'support_count':SUP[n],'support_probability':str(Fraction(SUP[n],A2[n]))})
z=json.loads((B/'independent-peaks-n6-r004.json').read_text());N=36
hist={int(k):int(v) for k,v in z['peak_histogram'].items()};dist={int(k):int(v) for k,v in z['distance_counts'].items()}
assert sum(hist.values())==factorial(N) and sum(dist.values())==A2[6]
assert hist[1]==int(z['A1']) and hist[2]==A2[6]
roots=list(map(int,z['single_peak_root_counts']));assert sum(roots)==hist[1]
for v in range(N):assert all(roots[v]==roots[transform(6,v,k)] for k in range(8))
Q=Fraction(int(z['squared_euclidean_total']),N*A2[6]);assert 0<Q<2
print(json.dumps({'verified':True,'beta_support':rows,'Q6':str(Q),'R6':str(Fraction(A2[6],hist[1]))},indent=2))
