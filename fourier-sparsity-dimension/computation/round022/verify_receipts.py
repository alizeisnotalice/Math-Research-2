from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent
read=lambda f:json.loads((p/f).read_text())
sha=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest()
for name in ['type2_core','gauss_and_cross_audit','cross_partition_audit','mm_cross_witness']:
 assert read(name+'.json')['source_sha256']==sha(name+'.py')
base=read('type2_core.json');assert base['sign_patterns']==16384 and base['allowed_count']==128 and base['orbit_representatives']==[22]
for host in ['local','macbook']:
 out='metal-'+host+'-result.json';d=read(out);r=read(host+'-receipt.json')
 assert r['source_sha256']==sha('R22Type2CoreMetal.swift') and r['output_sha256']==sha(out)
 assert r['compile_exit_code']==r['run_exit_code']==0
 assert d['legal_sign_masks_14bit']==base['allowed_masks']
 assert d['signed_walsh_coefficients_checked']==1048576 and d['candidate_sign_masks_checked']==16384
 assert d['gpu_cpu_mismatch_count']==0 and d['all_checks_passed']
 assert d['invalid_cpu_coefficients']==d['invalid_gpu_coefficients']==782336
assert all(x['both_bent'] for x in read('mm_cross_witness.json')['rows'])
print('PASS: source bindings, both GPU results, complete mask agreement, exact MM witnesses')
