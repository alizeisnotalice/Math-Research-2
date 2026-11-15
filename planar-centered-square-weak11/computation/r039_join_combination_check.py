#!/usr/bin/env python3
"""Bind independently audited R39 enumeration and sweep for one saved segment."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sweep_path,enum_path,review_path,note_path,out=map(lambda x:Path(x).resolve(),sys.argv[1:])
sweep=json.loads(sweep_path.read_text());enum=json.loads(enum_path.read_text());note=json.loads(note_path.read_text());text=review_path.read_text()
r=sweep['real'];assert sweep['status']==enum['status']==note['status']=='passed'
assert enum['slab_sha256']==note['slab_sha256']==r['slab_sha256']
assert note['slope_sha256']==r['slope_sha256']
assert enum['module_sha256']==r['slab_code_sha256']
assert enum['all_active_box_masses_radii_endpoints_verified']==r['boxes']
assert F(enum['H'])==F(r['H']) and F(enum['delta_empty'])==F(r['delta_empty'])
assert enum['slab_sha256'] in text and enum['module_sha256'] in text and '结论：PASS' in text
assert F(r['delta'])==min(F(r['H']),F(r['delta_empty']),F(r['endpoint_radius']))>0
assert 0<F(r['beta_lower'])<=F(r['beta_upper'])
assert not enum['all_rejection_classes_individually_replayed']
result={'status':'passed_combined_local_positive_segment',
 'scope':'enumerator code/completeness review plus independent real active-endpoint and y-union replay; not independent replay of all 44M rejected boxes',
 'raw_join_conditional_scope_preserved':True,'slab_sha256':r['slab_sha256'],'slope_sha256':r['slope_sha256'],
 'claim':'A(s)=A(s0)+beta*(s-s0) for |s-s0|<=delta, with the recorded strictly positive beta interval',
 'beta_lower':r['beta_lower'],'beta_upper':r['beta_upper'],'delta':r['delta'],
 'normalization':r['normalization'],'limiter_is_not_a_certified_area_saturation_wall':True,
 'evidence_sha256':{str(p):sha(p) for p in [sweep_path,enum_path,review_path,note_path]},
 'script_sha256':sha(Path(__file__)),'actual_command':[sys.executable,str(Path(__file__).resolve()),*map(str,[sweep_path,enum_path,review_path,note_path,out])]}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'script_sha256':result['script_sha256'],'output_sha256':sha(out)}))
