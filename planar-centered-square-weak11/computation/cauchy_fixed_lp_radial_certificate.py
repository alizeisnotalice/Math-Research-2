"""Exact certificate for radial monotonicity, not for the boundary minimum."""
import hashlib,json,pathlib,sys
from fractions import Fraction as Q
base=pathlib.Path(__file__).resolve().parent
p=base/'inputs/cauchy-smoothed-interface/fixed-lp-input.json'
x=json.loads(p.read_text())
A=Q();B=Q();rows=[]
for i,o in enumerate(x['orbits']):
    c=Q(o['coefficient']);R=Q(o['R'])
    assert c>=0 and R>1
    points=[tuple(Q(a) for a in z) for z in o['points']]
    assert len(set(points))==len(points)
    assert all(sum(map(abs,z))==R for z in points)
    assert set(points)=={tuple(-a for a in z) for z in points}
    A+=c;B+=c/(R-1)**2
    rows.append(dict(orbit=i,R=str(R),point_count=len(points),coefficient_exact_binary64=str(c),radial_B_contribution=str(c/(R-1)**2)))
margin=A-B
assert margin>0
out=dict(python=sys.version,input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),interpretation='each nonnegative coefficient multiplies U minus uniform orbit average V_t; binary64 JSON values interpreted as exact rationals',scope='strict full-D radial monotonicity for all common Cauchy smoothing times t>=0; NOT positive boundary minimum or mass/integral certificate',A_exact=str(A),A_float=float(A),B_exact=str(B),B_float=float(B),margin_exact=str(margin),margin_float=float(margin),normalized_B_over_A_exact=str(B/A),normalized_B_over_A_float=float(B/A),metadata_total_mass_exact=str(Q(x['total_mass'])),A_minus_metadata_exact=str(A-Q(x['total_mass'])),rows=rows)
print(json.dumps(out,indent=2))
