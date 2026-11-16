#!/usr/bin/env python3
"""Integer consistency checks for all fixed-point mean/variance ledgers."""
import json
from fractions import Fraction
from pathlib import Path

expected={3:(143808,633600),4:(779299174080,5272309842240),5:(7063261481362958280000,67631454448823894520000)}
rows=[]
for n,(W,D) in expected.items():
    data=json.loads(Path(f'mean-variance-n{n}.json').read_text())
    t={k:int(v) for k,v in data['totals'].items()}
    assert (t['W'],t['D'])==(W,D)
    for prefix in ('mean_square','within','between'):
        assert t[prefix+'_lower_numerator']<=t[prefix+'_upper_numerator']
    ledger=data['rank_minbranch']
    assert sum(int(x['W']) for x in ledger)==W
    assert sum(int(x['D']) for x in ledger)==D
    den=n*n*W*(int(data['fixed_point_L'])**2)
    wi=Fraction(t['within_lower_numerator'],den),Fraction(t['within_upper_numerator'],den)
    be=Fraction(t['between_lower_numerator'],den),Fraction(t['between_upper_numerator'],den)
    q=Fraction(D,n*n*W)
    assert wi[0]+be[0]<=q<=wi[1]+be[1]
    rows.append({'n':n,'within_interval':[str(x) for x in wi],'between_interval':[str(x) for x in be]})
out={'status':'passed','rows':rows}
Path('verify-mean-variance.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
