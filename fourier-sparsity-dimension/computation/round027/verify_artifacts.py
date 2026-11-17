from pathlib import Path
import json,hashlib
p=Path(__file__).parent
for src,out in [('check_fourier_sparsity_le4.py','fourier_sparsity_le4_result.json'),('check_twin_closure_and_pairings.py','twin_closure_and_pairings_result.json')]:
 d=json.loads((p/out).read_text())
 assert d['script_sha256']==hashlib.sha256((p/src).read_bytes()).hexdigest()
a=json.loads((p/'fourier_sparsity_le4_result.json').read_text())
assert a['enumerated_functions']==65536 and a['failure_count']==0
assert a['support_histogram']=={'1':32,'4':1120,'8':3840,'10':26880,'16':33664}
a=json.loads((p/'twin_closure_and_pairings_result.json').read_text())
for d in a['fixed_a_one_pair_codes']:
 q=1<<d['k'];m=q//2
 assert d['D_dim']==d['D4_rank']==m-1 and d['D_size_enumerated']==1<<(m-1)
 assert d['D4_weight4_count']==m*(m-1)//2 and d['H_dim']==m+1
 assert d['RM1_in_H'] and d['hypergraph_components']==1
 assert d['twin_classes']==[[2*i,2*i+1] for i in range(m)]
assert a['rm1_k3_twin_classes']==[[i] for i in range(8)]
e,o=a['k4_models']
assert e['twin_classes'][:3]==[[0,1],[2,3],[4,5]]
assert o['twin_classes'][:3]==[[0],[1,2,3],[4,8,12]]
assert all(d['RM1_in_H'] and d['D_dim']==d['D4_rank']==2 for d in (e,o))
d=a['k3_pair_block_walsh']
assert d['bent_F_count']==896 and d['cases_checked']==8064 and d['blocks_per_case']==4
assert d['all_passed'] and d['first_failure'] is None and d['coverage_histogram']=={'1':516096}
print('PASS: final source bindings, finite twin models and 8064 pair-block spectra')
