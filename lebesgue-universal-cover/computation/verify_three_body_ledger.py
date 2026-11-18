"""Verify the two stages' coverage counts and the rational demo weighting."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

base = Path(__file__).resolve().parent
stage = (base/'040-refinement.json').read_bytes()
first = json.loads(stage)
second = json.loads((base/'040-closure.json').read_text())
example = json.loads((base/'040-dual-example.json').read_text())
target = F(3339,4000)
assert F(first['target']) == F(second['target']) == target
assert second['source_stage_sha256'] == hashlib.sha256(stage).hexdigest()
expected = set()
for p in first['parts']:
    c = p['counts']
    terminals = len(p['unresolved']) + c.get('taylor',0) + c.get('corners',0)
    assert c.get('checked',0) == 2*terminals-c.get('seed_tiles',0)
    expected.update((p['part'],i) for i in range(len(p['unresolved'])))
covered = set()
for entry in second['proofs']:
    pair = entry['part'],entry['residual_index']
    assert pair not in covered
    covered.add(pair)
    proof = entry['proof']
    assert proof['method'] in ('single','dual') and F(proof['bound'])>target
    if proof['method']=='dual':
        assert 0<=F(proof['weight_first'])<=1
assert covered==expected and len(covered)==first['unresolved']
a,b=example['corner_lower_integer_rows']
w,v=map(F,example['weights'])
assert len(a)==len(b)==16 and min(w,v)>=0 and w+v==1
assert min(w*x+v*y for x,y in zip(a,b))/example['scale']==F(example['combined_lower'])>target
print('Verified full terminal counts, all remaining boxes, and the 16 rational demo inequalities.')
