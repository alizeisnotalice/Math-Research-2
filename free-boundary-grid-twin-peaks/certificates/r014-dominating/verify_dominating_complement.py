#!/usr/bin/env python3
"""Integer checks for the n=3,4,5 dominating-complement ledgers."""
import json
from pathlib import Path

expected={3:(143808,633600),4:(779299174080,5272309842240),5:(7063261481362958280000,67631454448823894520000)}
rows=[]
for n,(A2,Dall) in expected.items():
    d=json.loads(Path(f'dominating-n{n}.json').read_text())
    cats=d['categories'];u=d['covered_union']
    assert int(u['W'])==sum(int(x['W']) for x in cats.values())
    assert int(u['D'])==sum(int(x['D']) for x in cats.values())
    assert int(u['W'])<=A2 and int(u['D'])<=Dall
    assert all(x['h_factorial_failures']==0 for x in cats.values())
    assert all(x['state_count']%x['c0']==0 and x['h_factorial_failures']==0 for x in d['partition_layers'])
    rr=d['selected_side_ratio']
    assert int(rr['less_than_one']['W'])+int(rr['at_least_one']['W'])==int(u['W'])
    assert rr['less_than_one']['events']+rr['at_least_one']['events']==u['events']
    rows.append({'n':n,'covered_W':u['W'],'covered_D':u['D'],'partition_layers':len(d['partition_layers'])})
out={'status':'passed','rows':rows}
Path('verify-dominating-complement.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
