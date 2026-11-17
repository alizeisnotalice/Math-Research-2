from pathlib import Path
import hashlib,json
p=Path(__file__).parent
d=json.loads((p/'d4_residual_bound_result.json').read_text())
assert d['script_sha256']==hashlib.sha256((p/'check_d4_residual_bound.py').read_bytes()).hexdigest()
for a,t in zip(d['sharp_cubic_models'],[4,8]):
 assert a['d']==a['row_rank']==t and a['edges']==3*t//2
 assert a['row_weights']==[4]*t and a['row_support_label_xors']==[0]*t
 assert a['matrix_column_count_N']==a['r']==5*t//2
 assert all(a[k] for k in ['all_matrix_columns_distinct','all_physical_point_labels_distinct','RM1_orthogonality_verified','all_twin_classes_singleton'])
a=d['random_affine_weight4_trials']
assert a['trials']==len(a['outcomes'])==80 and a['seed']==20261117 and a['q']==32
assert a['all_N_bounds_pass'] and a['all_r_bounds_pass']
for t in a['outcomes']:
 assert t['D_rank']==t['D4_rank'] and t['D_size']==1<<t['D_rank']
 assert t['N']<=4*t['D4_rank'] and 2*t['r']<=5*t['D4_rank']
 assert t['N']+t['r']<=5*t['D4_rank']
print('PASS: final source binding, two sharp embeddings and 80 finite relation models')
