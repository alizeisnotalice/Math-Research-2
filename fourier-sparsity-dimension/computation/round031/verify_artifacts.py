from pathlib import Path
import hashlib,json
p=Path(__file__).parent
d=json.loads((p/'r31_scalar_rank4_absorption_result.json').read_text())
assert d['script_sha256']==hashlib.sha256((p/'check_r31_scalar_rank4_absorption.py').read_bytes()).hexdigest()
for k,v in {'q':16,'base_rank4_graph_blocks':4,'code_rank':9,'phase_count':512,'residual_size':128,'M':8,'block_size':32,'gamma_tiles':8,'gamma_tile_size':16,'fixed_D_size':128,'joint_labelled_translation_rank':5,'block_constant_kernel_rank':4,'freeze_rank_mod_scalar_Q':0,'tau_x5_parameter_mask':64}.items():assert d[k]==v
for k in ['all_phases_bent','fixed_D_is_u_frequency_zero','partial_gamma_tiles_fixed_disjoint_pm_q','translation_labels_linear_all_parameters','translation_kernel_equals_block_constants','per_block_quadratic_variation_scalar_Q']:assert d[k]
print('PASS: final source binding, 512 phases and exact scalar-Q absorption with fixed full patch partition')
