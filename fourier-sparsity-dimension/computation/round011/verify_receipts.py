#!/usr/bin/env python3
"""Bind current sources, historical seed inputs, two-machine and independent audits."""
import argparse,hashlib,json,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
base=Path(__file__).parent;root=base.parents[2]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
a=json.loads((base/'period-local-final.json').read_text())
b=json.loads((base/'period-remote-final.json').read_text())
g=json.loads((base/'graph-spectra-verification.json').read_text())
assert a['status']==b['status']==g['status']=='passed'
assert a['source_sha256']==b['source_sha256']==sha(base/'audit_common_affine_derivative_periods.py')
assert g['source_sha256']==sha(base/'verify_graph_spectra.py')
assert a['audit_input_sha256']==b['audit_input_sha256']==g['input_sha256']==sha(base/'period-input.json')
assert a['rows']==b['rows'] and len(a['rows'])==len(g['rows'])==20
assert sum(x['P_dimension']==0 for x in a['rows'])==10
for x,y in zip(a['rows'],g['rows']):
    assert x['row_index']==y['row_index'] and x['source_sha256']==y['source_sha256']
    assert x['P_dimension']==y['period_dimension']
    assert sorted(x['P_elements'])==sorted((v>>1)&63 for v in y['periods'])
    if x['P_dimension']==0:assert x['all_32768_slices_independently_bent']
    relative='fourier-sparsity-dimension/'+x['source_file'].split('fourier-sparsity-dimension/',1)[1]
    assert sha(root/relative)==x['source_sha256']
c=json.loads((base/'rank-and-code-bound-controls.json').read_text())
d=json.loads((base/'dim7-model-and-row-obstruction.json').read_text())
assert c['status']=='passed' and c['source_sha256']==sha(base/'verify_rank_and_code_bounds.py')
assert d['runtime']['script_sha256']==sha(base/'find_dim7_local_model.py')
assert d['single_boolean_row_compatible_phase_count']==0
out={'status':'passed','period_records_equal_across_machines':20,'independent_output_domain_checks':20,
     'aperiodic_examples':10,'source_and_seed_hashes_verified':True,
     'local_slice_seconds':a['elapsed_seconds'],'remote_slice_seconds':b['elapsed_seconds'],
     'direct_graph_seconds':g['elapsed_seconds'],'python':sys.version,
     'source_sha256':sha(Path(__file__)),
     'input_sha256':{n:sha(base/n) for n in ['period-input.json','period-local-final.json','period-remote-final.json','graph-spectra-verification.json','rank-and-code-bound-controls.json','dim7-model-and-row-obstruction.json']}}
args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
