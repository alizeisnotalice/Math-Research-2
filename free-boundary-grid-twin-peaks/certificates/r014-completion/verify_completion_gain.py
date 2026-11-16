#!/usr/bin/env python3
"""Consistency checks for the exact n=3,4 completion-gain ledgers."""
import json
from fractions import Fraction
from pathlib import Path
expected={3:(143808,633600),4:(779299174080,5272309842240)}
rows=[]
for n,(A2,Dall) in expected.items():
    d=json.loads(Path(f'completion-n{n}.json').read_text());u=d['available_union'];av=d['availability']
    assert int(u['W'])==sum(int(x['W']) for x in av.values())
    assert int(u['D'])==sum(int(x['D']) for x in av.values())
    assert int(d['all_pairs']['W'])==int(u['W'])+int(d['unavailable_both_H_zero']['W'])
    assert int(d['all_pairs']['D'])==int(u['D'])+int(d['unavailable_both_H_zero']['D'])
    assert int(d['all_pairs']['W']) == A2 and int(d['all_pairs']['D']) == Dall
    assert int(u['W'])<=A2 and int(u['D'])<=Dall
    gp=d['gamma_union_positive_exponents'];gn=d['gamma_union_negative_exponents']
    assert int(gn['le_4_neg_g']['W'])<=int(gn['le_2_neg_g']['W'])<=int(gp['le_1']['W'])<=int(gp['le_2_pow_g']['W'])<=int(gp['le_4_pow_g']['W'])<=int(u['W'])
    assert sum(int(x['W']) for x in d['rank_minbranch'])==int(u['W'])
    assert sum(int(x['D']) for x in d['rank_minbranch'])==int(u['D'])
    q=Fraction(int(u['D']),n*n*int(u['W']))
    iv=d['sum_F_over_Dstar_rigorous_interval'];assert int(iv['lower_numerator'])<int(iv['upper_numerator'])
    assert int(iv['lower_numerator']) > int(iv['denominator'])
    assert int(iv['upper_numerator'])-int(iv['lower_numerator']) == u['pairs']
    rows.append({'n':n,'coverage':str(Fraction(int(u['W']),A2)),'conditional_Q':str(q)})
out={'status':'passed','rows':rows}
Path('verify-completion-gain.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
