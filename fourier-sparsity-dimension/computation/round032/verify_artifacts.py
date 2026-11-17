"""Validate bindings and independently check the small projection certificates.

This does not rerun the published affine classification or the full kernel counts.
"""
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
read=lambda n:json.loads((P/n).read_text())
sha=lambda n:hashlib.sha256((P/n).read_bytes()).hexdigest()
r=read('rm48_weight32_kernel_audit.json');reps=read('rm48_weight32_representatives.json')['representatives']
assert r['script_sha256']==sha('audit_rm48_kernels.py')
assert r['input_sha256']==sha('rm48_weight32_representatives.json')
counts=[32768,8192,2048,3072,512,0,1920,1344,192]
for rep,row,count in zip(reps,r['representatives'],counts,strict=True):
 assert row['id']==rep['id'] and row['avoiding_kernel_count']==count
 assert row['kernel_subspaces_checked']==97155
 masks=[sum(1<<i for i in mon) for mon in rep['monomials_zero_based']]
 A={x for x in range(1<<rep['affine_span_dimension']) if (rep['constant_term']+sum((x&v)==v for v in masks))%2}
 assert len(A)==32
 D={x^y for x in A for y in A}
 assert len(D)==row['difference_size']
 if count:
  K={0}
  for v in row['first_rref_basis']:K|={x^v for x in list(K)}
  assert len(K)==8 and K&D=={0}
 else:
  assert rep['id']=='span7_code11' and D==set(range(128))
  # Every 3-space in F2^8 intersects this 7-space in dimension at least 2.
a=read('transvection_avoidance_result.json')
assert a['script_sha256']==sha('check_transvection_avoidance.py')
g=a['k3_geometry'];assert g['equivalence_checks_passed']==12870 and g['degree_le_2_B_count']==870 and g['degree_le_2_all_projectable']
assert g['projectable_B']==2790 and g['nonprojectable_B']==10080
assert all(a['k4_actual_bent_transvection'].values())
for source,result,inp,minimum in [('k5_transvection_candidate_count.py','k5-transvection-candidates.json','k5-transvection-input.json',192),('r18_tensor10_transvection_histogram.py','r18-tensor10-transvection-histogram.json','r18-tensor10-transvection-input.json',128)]:
 d=read(result)
 assert d['source_sha256']==sha(source) and d['input_sha256']==sha(inp)
 assert d['transvections_tested']==522753 and sum(d['weight_histogram'].values())==522753
 assert d['weight64_candidate_count']==0 and d['minimum_positive_weight']==minimum
print('PASS: source/input bindings, eight explicit kernels, full-hyperplane obstruction, geometry and two seed histograms')

# Both complete sign scans use the same eight physical gauge points but
# different support orderings; compare the actual negative point sets.
t=read('exceptional_signed_mask_result.json');c=read('signed_plateaued_enum_result.json')
assert t['source_sha256']==sha('check_exceptional_signed_mask.py')
assert c['source_sha256']==sha('signed_plateaued_enum.c')
assert t['checked']==c['states_checked']==1470943
assert len(t['solutions'])==c['solution_count']==1
neg={x for i,x in enumerate(t['points']) if t['solutions'][0]>>i&1}
assert neg=={x for i,x in enumerate(c['support']) if c['negative_sign_mask']>>i&1}
W=[sum((-1 if x in neg else 1)*(-1 if (a&x).bit_count()%2 else 1) for x in t['points']) for a in range(128)]
assert set(W)=={-16,0,16} and sum(v!=0 for v in W)==16
print('PASS: independent complete sign scans agree; witness has 16 nonzero Walsh coefficients of magnitude 16')
