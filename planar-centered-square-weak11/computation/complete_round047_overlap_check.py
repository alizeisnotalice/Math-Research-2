#!/usr/bin/env python3
"""Bind R47 gap/slope replay to candidate/outer-Y/upstream-domain acceptance."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
own_path,support_path,review_path,old_path,out=map(lambda s:Path(s).resolve(),sys.argv[1:])
own=json.loads(own_path.read_text());support=json.loads(support_path.read_text());old=json.loads(old_path.read_text())
assert own['status']=='passed_gap_and_slope_replay_pending_catalogue_binding' and support['status']=='passed'
assert own['contract_sha256']==support['math_sha256']=='3e391beace80ec39a668eaeb81de218d8ab2d982da39d919f8f3ead749d5a3b7'
assert own['module_sha256']==support['module_sha256']=='929cf92baca532b92159ff81542390acbff3688386dd5104decd0ee55a8ed524'
assert own['index_sha256']==support['six_index_sha256'] and own['new_center_index_sha256']==support['five_new_catalogue_index_sha256']
assert sha(review_path)=='070934c47ef6c5f676f1e2b5aad4f270bb4ae7e425fff77f7c4bbca6b5036f1c'
assert sha(old_path)=='e9732e787d4ae16f1a087f2a3f9274449f2acdbc0baa94b6eeb43de6762b2504'
a={r['name']:r for r in own['real_radius_replays']};b={r['name']:r for r in support['six_support_attachments']}
assert set(a)==set(b) and len(a)==6
for name,r in a.items():
    s=b[name]
    for key in ['source_slab_sha256','output_sha256']:assert r[key]==s[key]
    assert r['boxes']==s['candidate_boxes'] and r['x_classes']==s['canonical_x_classes']
    assert r['positive_outer_y_strips']==s['positive_y_strips']
    assert s['Fx_remains_active_and_references_balance'] and s['limiter_owner_strip_verified']
    assert F(r['final_delta'])==F(s['final_delta'])==min(F(r['H']),F(r['delta_empty']),F(r['overlap_radius']))
    assert F(r['improvement_factor'])==F(s['improvement_factor'])>=1
    assert F(r['beta_lower'])>0 or F(r['beta_upper'])<0
new={r['direction']:r for r in own['new_center_slope_replays']};cats={r['name']:r for r in support['five_new_catalogues']}
assert set(new)==set(cats)==set(a)-{'r039'}
for name,r in new.items():
    c=cats[name];assert r['hashes']['slab.json']==c['slab_sha256']==a[name]['source_slab_sha256']
    assert r['boxes']==c['active_boxes_verified'] and r['velocity_nonzero_labels']==c['labels']
    for x,y in [('s0_physical','physical_s0'),('s0_normalized','normalized_s0'),('H','H'),('delta_empty','delta_empty')]:assert F(r[x])==F(c[y])
assert sum(r['boxes'] for r in new.values())==support['total_new_active_boxes_verified']==133259
assert old['status']=='passed_combined_local_positive_segment' and old['slab_sha256']==a['r039']['source_slab_sha256']
assert F(old['beta_lower'])==F(a['r039']['beta_lower']) and F(old['beta_upper'])==F(a['r039']['beta_upper'])
assert {name for name,r in a.items() if F(r['improvement_factor'])>1}=={'r039','q0'}
assert not support['large_bbox_loops_replayed'] and not support['recenter_remaining_domain_automatically_validated']
result={'status':'passed_combined_R47_six_radii_and_five_new_slopes',
 'scope':'six original-slab candidate families with independently accepted supports/domains; five new center slopes genuinely replayed, R39 slope reused',
 'contract_sha256':own['contract_sha256'],'module_sha256':own['module_sha256'],
 'six_index_sha256':own['index_sha256'],'five_new_index_sha256':own['new_center_index_sha256'],
 'claim':'Each accepted complete union area remains affine with its recorded beta on |s-s0|<=final_delta; only r039 and q0 final radii improve.',
 'radii_and_slopes':list(a.values()),'recenter_external_remaining_domain_premise_preserved':True,
 'nearest_effective_event_or_saturation_claimed':False,'R39_R45_frozen_outputs_modified':False,
 'evidence_sha256':{str(p):sha(p) for p in [own_path,support_path,review_path,old_path]},
 'script_sha256':sha(Path(__file__)),'actual_command':[sys.executable,str(Path(__file__).resolve()),*map(str,[own_path,support_path,review_path,old_path,out])]}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'script_sha256':result['script_sha256'],'output_sha256':sha(out)}))
