from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent
read=lambda f:json.loads((p/f).read_text())
for name in ['four_relation_closure','nonaffine_component_spectrum']:
 assert read(name+'.json')['source_sha256']==hashlib.sha256((p/(name+'.py')).read_bytes()).hexdigest()
a=read('four_relation_closure.json');assert a['small_seed_count']==67 and a['all_pair_partitions_equal']
assert len(a['small_proper_closures'])==1
assert [(r['seed_dimension'],r['closure_dimension']) for r in a['large_models']]==[(14,15),(14,14),(15,16)]
assert sorted(map(len,a['large_models'][1]['components']))==[1]*10+[6]
b=read('nonaffine_component_spectrum.json');assert b['all_H_phases_checked']==128 and b['all_bent']
assert b['frequency_component']==[0,1,2,4] and not b['frequency_component_affine'] and b['partial_support_size']==32
print('PASS: source bindings, all 67 seeds and actual nonaffine partial spectrum')
