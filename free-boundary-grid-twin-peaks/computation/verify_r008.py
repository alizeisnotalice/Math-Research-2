#!/usr/bin/env python3
"""Audit archived exact ledgers; does not replace enumeration/proofs."""
from pathlib import Path
from fractions import Fraction
import json
ROOT=Path(__file__).resolve().parent
cpp=json.loads((ROOT/'results/n8-kernel-counterexample-cpp-r008.json').read_text())
py=json.loads((ROOT/'results/n8-kernel-counterexample-python-r008.json').read_text())
a=py['n8']
h=int(cpp['hS']); H=int(cpp['p_numerator'])
assert (h,H)==(a['hS'],a['H'])
assert a['pair_total']==49*48*H
assert a['independent_sets']==cpp['independent_sets']==3843618
assert sum(x['independent_sets'] for x in a['independent_layers'])==3843618
assert 14*H-4*h==py['comparison_gap_boundary_H_minus_4hS']>0
q=Fraction(14*H,h)
assert q==Fraction(198947127410665164977798511882264205133,49059448277247751015473936500562374772)
amp=Fraction(34*49*48,69*68)*Fraction(H,h)
assert amp==Fraction(198947127410665164977798511882264205133,40298832513453509762710733554033379277)
assert amp>q>4
# Verify geometry directly, independently of integer counting.
S={tuple(x) for x in cpp['S_cells']}
def adj(v):
    x,y=v
    return {(a,b) for a,b in [(x-1,y),(x+1,y),(x,y-1),(x,y+1)] if 0<=a<8 and 0<=b<8}
deg={v:len(adj(v)&S) for v in S}
assert sorted(deg.values())==[1,1]+[2]*13
assert {v for v,d in deg.items() if d==1}=={(7,0),(0,7)}
seen={next(iter(S))}
while True:
    new=seen|set().union(*(adj(v)&S for v in seen))
    if new==seen:break
    seen=new
assert seen==S
assert len(set().union(*(adj(v) for v in S))-S)==14
manifest=json.loads((ROOT/'certificates/r008-kernel-selected/n7-critical-d4/compact-manifest.json').read_text())
assert len(manifest['runs'])==30 and all(x['exit_code']==0 and x['candidate_holds'] for x in manifest['runs'])
print(json.dumps({'status':'pass','independent_integer_certificates_agree':True,'gap':14*H-4*h,'n8_scaled_rational':str(q),'amplified_rational':str(amp),'geodesic':True,'n7_orbits':30},indent=2))
