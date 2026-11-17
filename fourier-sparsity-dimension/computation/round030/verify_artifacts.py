from pathlib import Path
import json,hashlib
p=Path(__file__).parent
d=json.loads((p/'r30_small_examples_result.json').read_text())
assert d['script_sha256']==hashlib.sha256((p/'check_r30_small_examples.py').read_bytes()).hexdigest()
assert d['q']==8
for key,dim,phases,size,trank,orank in [('four_tile',7,128,32,3,4),('three_tile',4,16,24,1,3)]:
 a=d[key]
 assert (a['dimension'],a['phase_count'],a['residual_size'],a['fixed_total_frequency_size'],a['translation_rank'],a['offset_rank'])==(dim,phases,size,size,trank,orank)
 assert all(a[k] for k in ['all_bent','blocks_size8','nonparallel','partial_tiles_fixed_disjoint_pm8','translation_kernel_offset_span'])
assert d['four_tile']['fixed_total_equals_u_frequency_zero_D'] and not d['four_tile']['fixed_total_nonaffine']
assert d['three_tile']['fixed_total_nonaffine']
print('PASS: final source binding, 144 bent phases, fixed nonparallel tiles and translation/kernel ranks')
