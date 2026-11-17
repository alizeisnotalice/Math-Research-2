"""Check archived source bindings and independent finite result agreement."""
from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parent
read=lambda f:json.loads((p/f).read_text())
sha=lambda f:hashlib.sha256((p/f).read_bytes()).hexdigest()
for src,out in [('six_point_signs.py','six-point-signs.json'),('six_point_fractional.py','six-point-fractional.json'),('enumerate_quadric_core_free1.py','quadric_core_free1.json'),('three_half_blocks_joint_audit.py','three-half-block-joint-audit.json')]:
 assert read(out)['source_sha256']==sha(src)
f=read('six-point-fractional.json'); expected={tuple(x) for x in f['allowed_numerators']}
assert f['grid_points']==531441 and len(expected)==33 and f['nonsaturated_count']==0
for host,source in [('local','R21LocalMetal.swift'),('macbook','R21MacBookMetal.swift')]:
 out='metal-'+host+'-result.json'; d=read(out); r=read(host+'-receipt.json')
 assert r['source_sha256']==sha(source) and r['output_sha256']==sha(out)
 assert r['compile_exit_code']==r['run_exit_code']==0
 assert d['coefficient_vectors_checked']==531441 and d['signed_walsh_coefficients_checked']==8503056
 assert d['gpu_cpu_mismatch_count']==0 and d['solution_count']==33
 assert d['invalid_cpu_coefficients']==d['invalid_gpu_coefficients']==7939296
 assert {tuple(x['m']) for x in d['solutions']}==expected
 assert d['all_checks_passed']
a=read('six-point-signs.json'); assert a['sign_patterns_tested']==64 and a['allowed_count']==32 and a['single_affine_orbit']
a=read('quadric_core_free1.json'); assert a['sign_tables']==4096 and a['local_plateaued_tables']==a['common_free_character_tables']==64
a=read('three-half-block-joint-audit.json'); assert a['phase_count']==128 and a['all_phases_bent'] and a['joint_cylinders_disjoint']
assert a['joint_cylinder_count']==7 and a['parameter_projection_rank']==a['normal_union_rank']==4
assert a['K_dimension']==3 and set(a['K_block_restriction_ranks'].values())=={1}
print('PASS: source bindings, both GPU coefficient checks, identical 33 solutions, free-variable and mixed-partition witnesses')
