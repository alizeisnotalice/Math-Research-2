from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent
read=lambda f:json.loads((p/f).read_text())
sha=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest()
for s,o in [('two_patch_audit.py','two_patch_audit.json'),('desarguesian_cross_partition_audit.py','desarguesian-cross-partition-audit.json')]:
 assert read(o)['source_sha256']==sha(s)
py=read('two_patch_audit.json');results=[]
for h in ['local','macbook']:
 out='metal-'+h+'-result.json';d=read(out);r=read(h+'-receipt.json')
 assert r['source_sha256']==sha('R23TwoPatchMetal.swift') and r['output_sha256']==sha(out)
 assert r['compile_exit_code']==r['run_exit_code']==0 and d['all_checks_passed']
 assert '-module-cache-path' in r['compile_command']
 for x,y in zip(py['rows'],d['rows']):
  for key in ['mu','free_dimension','phase_count','patch_sizes','joint_patch_codimensions','parameter_projection_rank','K_dimension','K_restriction_rank']:assert x[key]==y[key]
  for key in ['gpu_cpu_mismatch_count','full_cpu_bad_absq','full_gpu_bad_absq','partial_patch_cpu_bad_absq','partial_patch_gpu_bad_absq']:assert y[key]==0
 assert len(d['rows'])==3 and sum(2*x['phase_count']*x['N'] for x in d['rows'])==1081344
 results.append(d['rows'])
assert results[0]==results[1]
s=read('desarguesian-cross-partition-audit.json');assert [x['generator_rank'] for x in s['rows']]==[11,17,27]
assert all(x['blocks_partition_domain'] and x['base_phase_bent'] for x in s['rows'])
print('PASS: final source/output bindings, both GPU rows, Python invariants, spread ranks')
