"""Independent arithmetic and prior-method cross-checks for round 5."""
from pathlib import Path
from fractions import Fraction
from math import factorial, comb
import json
from verify_independent_peak_dp import EXPECTED
B=Path(__file__).parent/'results'
def read(name): return json.loads((B/name).read_text())
for n in range(3,7):
    z=read(f'd4-generic-n{n}-r005.json')
    a,b,d,s=EXPECTED[n]
    assert (int(z['A1']),int(z['A2']),int(z['squared_euclidean_total']))==(a,b,s)
    assert {int(k):int(v) for k,v in z['distance_counts'].items()}==d
    if n==6:
        old=read('independent-peaks-n6-r004.json')
        for key in ('A1','A2','peak_histogram','distance_counts','squared_euclidean_total'):
            assert old[key]==z[key]
z=read('n7-exact-r005.json'); N=49
hist={int(k):int(v) for k,v in z['peak_histogram'].items()}
assert all(v>0 for v in hist.values()) and sum(hist.values())==factorial(N)
assert sum(map(int,z['distance_counts'].values()))==hist[2]==int(z['A2'])
assert hist[1]==int(z['A1'])
n=7
assert Fraction(sum(k*v for k,v in hist.items()),factorial(N))==Fraction(3*n*n+3*n+2,15)
assert Fraction(sum(comb(k,2)*v for k,v in hist.items()),factorial(N))==(Fraction(n**4,50)+Fraction(n**3,25)-Fraction(11*n*n,450)-Fraction(n,15)-Fraction(13,315))
burn=read('d4-burnside-r005.json')
assert burn['verified'] and burn['D4_orbits_by_size']==z['representatives_by_size']
assert sum(burn['raw_layers'])==z['raw_independent_sets']==1280128950
rect=read('rectangle-summary-r005.json'); assert int(rect['A2_7'])==hist[2]
L={int(k):int(v) for k,v in rect['L'].items()}
C={j:comb(49,7*j)*L[j]*L[7-j] for j in range(2,6)}
assert C=={int(k):int(v) for k,v in rect['C'].items()}
assert sum(C.values())==int(rect['S'])
assert Fraction(rect['parity_over_A2'])==Fraction(max(sum(v for j,v in C.items() if j%2==q) for q in (0,1)),hist[2])
beta=B.parent/'certificates/r005-beta'
historical=json.loads((beta/'beta-zero-layered-n7.json').read_text())
seen=set();total=0
for i,o in enumerate(historical['edge_orbits']):
    v,w=o['representative'];assert abs(v%7-w%7)+abs(v//7-w//7)==1
    def trans(v,k):
        x,y=v%7,v//7
        if k>=4:x=6-x
        for _ in range(k%4):x,y=6-y,x
        return 7*y+x
    orbit={tuple(sorted((trans(v,k),trans(w,k)))) for k in range(8)}
    assert len(orbit)==o['multiplicity'] and not seen.intersection(orbit)
    seen.update(orbit)
    replay=json.loads((beta/f'orbits/orbit-{i:02d}.json').read_text())['edge_orbits'][0]
    for key in ('representative','multiplicity','E_K2','v_high','beta0_one_direction'):
        assert replay[key]==o[key]
    b0=int(o['beta0_one_direction']);assert b0==int(o['E_K2'])-int(o['v_high'])>=0
    total+=2*len(orbit)*b0
assert len(seen)==84 and total==int(historical['beta0_total'])
assert hist[2]-total==5548507403524089636774340128876992769846373065420288
q=Fraction(int(z['squared_euclidean_total']),49*hist[2]);r=Fraction(hist[2],hist[1])
print(json.dumps({'verified':True,'n3_to_n6_prior_method_full_distance_checks':True,'n6_full_histogram_check':True,'n7_histogram_moments_and_burnside':True,'Q7':str(q),'R7':str(r),'rectangle_far_lower_bound':rect['parity_over_A2']},indent=2))
