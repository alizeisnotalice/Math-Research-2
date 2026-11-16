#!/usr/bin/env python3
"""Exact integer consistency checks for the compact r015 ledgers."""
import json
from fractions import Fraction
from pathlib import Path
expected={3:(143808,633600),4:(779299174080,5272309842240),5:(7063261481362958280000,67631454448823894520000)}
rows=[]
for n,(A2,D) in expected.items():
 d=json.loads(Path(f'n{n}-r015.json').read_text()); a=d['all_first_merge_states']; u=d['available_union']; av=d['availability']; rem=d['unavailable_both_H_zero_weighted_remainder']
 assert (int(a['W']),int(a['D']))==(A2,D)
 assert int(u['W'])==sum(int(x['W']) for x in av.values())
 assert int(u['D'])==sum(int(x['D']) for x in av.values())
 assert int(u['W'])+int(rem['W'])==A2 and int(u['D'])+int(rem['D'])==D
 assert sum(int(x['W']) for x in d['rank_minbranch'])==int(u['W'])
 assert sum(int(x['D']) for x in d['rank_minbranch'])==int(u['D'])
 gp=d['gamma_union_positive_exponents']; gn=d['gamma_union_negative_exponents']
 assert int(gn['le_4_neg_g']['W'])<=int(gn['le_2_neg_g']['W'])<=int(gp['le_1']['W'])<=int(gp['le_2_pow_g']['W'])<=int(gp['le_4_pow_g']['W'])==int(u['W'])
 iv=d['sum_F_over_Dstar_rigorous_interval']; lo=Fraction(int(iv['lower_numerator']),int(iv['denominator'])); hi=Fraction(int(iv['upper_numerator']),int(iv['denominator'])); assert lo<hi
 if n==5:
  # Human-readable strict enclosure; the exact dyadic interval remains in JSON.
  assert Fraction(1573345373149118,10**9)<lo<hi<Fraction(1573345373149119,10**9)
 rows.append({'n':n,'coverage':str(Fraction(int(u['W']),A2)),'Q_all':str(Fraction(D,n*n*A2)),'Q_available':str(Fraction(int(u['D']),n*n*int(u['W']))),'sum_lower_gt_one':int(iv['lower_numerator'])>int(iv['denominator'])})
out={'status':'passed','rows':rows}
Path('verify-r015.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
