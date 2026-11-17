from pathlib import Path
import hashlib,json
p=Path(__file__).parent
d=json.loads((p/'r29_counterexamples_and_permutations_result.json').read_text())
assert d['script_sha256']==hashlib.sha256((p/'check_r29_counterexamples_and_permutation_families.py').read_bytes()).hexdigest()
assert [a['k'] for a in d['two_block_examples']]==[2,3,4,5,6]
for a in d['two_block_examples']:
 assert a['rank']==3 and a['phases']==8 and a['all_bent'] and a['failure'] is None
 assert a['local_budget_sum']==2 and a['supports_change'] and a['R_size']==2*a['q']
for a in d['r4_parallel_examples']:
 assert a['rank']==7 and a['phases']==128 and a['all_bent'] and a['failure'] is None
 assert len(a['slope_images'])==8 and all(sorted(s)==[0,1,2,3] for s in a['slope_images'])
assert [a['k'] for a in d['r4_parallel_examples']]==[3,4]
a=d['small_permutation_direction_enumeration']
assert [x['S'] for x in a]==[[0,1,2],[0,1,2,3],[0,1,2,4],[0,1,2,4,8]]
assert [x['max_dimension'] for x in a]==[1,3,2,2]
assert [x['legal_linear_subspace_count'] for x in a]==[6,60,27,155]
assert [x['permutation_count'] for x in a]==[6,24,24,120]
print('PASS: source binding, 296 bent phases and four exhaustive finite permutation-family checks')
