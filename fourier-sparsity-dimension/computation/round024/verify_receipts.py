from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent
read=lambda f:json.loads((p/f).read_text())
sha=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest()
for name in ['expurgate_small_directions','explicit_two_checks','degree_capped_triangle']:
 assert read(name+'.json')['source_sha256']==sha(name+'.py')
for host in ['local','macbook']:
 out='metal-'+host+'-result.json';d=read(out);r=read(host+'-receipt.json')
 assert r['source_sha256']==sha('R24ExpurgateMetal.swift') and r['output_sha256']==sha(out)
 assert r['input_sha256']==d['input_sha256']==sha('expurgate_small_directions.json')
 assert r['compile_exit_code']==r['run_exit_code']==0 and d['all_checks_passed']
 assert d['basis_rank']==15 and d['check_rank']==2 and d['retained_phase_count']==8192
 assert d['signed_walsh_coefficients_checked']==524288 and d['minimum_nonzero_direction_weight']==16
 assert d['gpu_cpu_mismatch_count']==d['cpu_bad_abs8']==d['gpu_bad_abs8']==0
 assert d['check_rows']==[15,224] and d['affine_basis_exact']
a=read('degree_capped_triangle.json');row=next(x for x in a['rows'] if (x['k'],x['ell'])==(4,2))
assert row['dimension']==21 and row['full_weight_enumeration_count']==2097152
hist={int(k):v for k,v in row['weight_histogram'].items()};assert sum(hist.values())==2097152 and min(k for k in hist if k)>0 and min(k for k in hist if k)==64
print('PASS: source/input/output bindings, both GPU certificates, complete 21-dimensional weight scan')
