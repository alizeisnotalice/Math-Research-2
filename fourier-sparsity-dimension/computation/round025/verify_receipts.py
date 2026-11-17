from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent
read=lambda f:json.loads((p/f).read_text())
sha=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest()
for name in ['sidon_seed_audit','weight_four_obstruction']:
 assert read(name+'.json')['source_sha256']==sha(name+'.py')
expected={int(x,16) for x in read('sidon_seed_audit.json')['rows'][0]['compatible_masks_hex']};assert len(expected)==384
for h in ['local','macbook']:
 out='metal-'+h+'-result.json';d=read(out);r=read(h+'-receipt.json')
 assert r['source_sha256']==sha('R25SidonMetal.swift') and r['output_sha256']==sha(out)
 assert r['compile_exit_code']==r['run_exit_code']==0 and d['all_checks_passed']
 assert d['all_truth_tables']==65536 and d['signed_walsh_coefficients_checked']==1048576
 assert d['gpu_cpu_mismatch_count']==0 and d['bent_count']==896
 assert {int(x,16) for x in d['compatible_masks_hex']}==expected and d['compatible_equals_MM']
 assert d['counterexample_affine_w_all_bent'] and d['counterexample_not_MM'] and d['counterexample_w_mask_0e_detects_illegal']
assert [x['compatible_phases_checked'] for x in read('weight_four_obstruction.json')['rows']]==[8,128,32768]
print('PASS: final source/output bindings, complete 384-mask equality and obstruction counts')
