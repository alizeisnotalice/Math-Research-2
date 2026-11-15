#!/usr/bin/env python3
"""R45 evidence binding and exact sign/bracket conclusions; no geometry run."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sweep_path,enum_path,enum_review_path,report_path,secant_path,prior_review_path,out=map(lambda x:Path(x).resolve(),sys.argv[1:])
sweep=json.loads(sweep_path.read_text());enum=json.loads(enum_path.read_text())
assert sweep['status']=='passed_sweep_replay_pending_catalogue_binding' and enum['status']=='passed'
assert sweep['index_sha256']==enum['index_sha256']=='855c330b1d36dc212956eefcfba840a831aeabe17d7f8ba6dfa94b5466986a88'
assert F(sweep['base_horizontal_period'])==F(enum['R8_base_period'])
left={r['name']:r for r in sweep['runs']};right={r['name']:r for r in enum['runs']}
assert set(left)==set(right) and len(left)==11
for name,r in left.items():
    e=right[name]
    assert r['hashes']==e['hashes']
    assert r['enum_sha256']==enum['source_module_sha256']
    assert r['input_sha256']==enum['R41_input_sha256']
    assert r['velocity_nonzero_labels']==e['labels']
    assert r['boxes']==e['active_boxes_verified']
    for a,b in [('s0_physical','physical_s0'),('s0_normalized','normalized_s0'),('H','H'),('delta_empty','delta_empty'),('tail_bound','tail_bound')]:assert F(r[a])==F(e[b])
    assert F(r['delta_physical'])==min(F(r['H']),F(r['delta_empty']),F(r['endpoint_radius']))>0
    assert (F(r['beta_lower'])>0 if r['sign']=='positive' else F(r['beta_upper'])<0)
assert sum(r['boxes'] for r in left.values())==enum['total_active_boxes_verified']==292849
assert not enum['all_rejected_boxes_individually_replayed']
assert sha(report_path)==enum['report_sha256']=='0017a70bee0bf505f37735150b363bf6fdb27d84f8ba0828f37c5b1677cc99d7'
assert sha(enum_review_path)=='7c95396316d0e3e41825c5d0b4b6a887b40a0e345464c8e36355553703596262'
assert sha(prior_review_path)=='13e9aaa0581e8333aebdccc0dac11e827b7f28b760014c3ff724b7fcb192ee23'
secant=json.loads(secant_path.read_text())
old=next(r for r in secant['directions']['(0, 1, 2, 3)']['segments'] if F(r['left'])==F(-11,14400000) and F(r['right'])==F(-7,9600000))
mid=left['negative-slope-all'];assert mid['sign']=='positive'
assert (F(old['left'])+F(old['right']))/2==F(mid['s0_normalized'])==F(-43,57600000)
assert F(old['slope_upper'])<0
assert len(sweep['five_direction_pairs'])==5
for pair in sweep['five_direction_pairs']:
    a,b=left[pair['left_run']],left[pair['right_run']]
    assert a['sign']=='positive' and b['sign']=='negative' and F(a['s0_normalized'])<F(b['s0_normalized'])
result={'status':'passed_combined_R45_eleven_local_segments',
 'scope':'eleven new R45 source/slab and sweep instances; source review and tests reused, no candidate-bbox classification loops rerun',
 'index_sha256':sweep['index_sha256'],'five_direction_pairs':sweep['five_direction_pairs'],
 'pair_conclusion':'each continuous one-dimensional direction has an interior maximum point or plateau between the recorded left-positive and right-negative centers',
 'no_uniqueness_exact_wall_or_multidimensional_claim':True,
 'extra_midpoint':{'run':mid['name'],'normalized_s0':mid['s0_normalized'],'beta_lower':mid['beta_lower'],'beta_upper':mid['beta_upper'],
   'prior_secant_left':old['left'],'prior_secant_right':old['right'],'prior_secant_upper':old['slope_upper'],
   'prior_secant_reused_not_recomputed':True,'not_monotone_decreasing_throughout_prior_secant_interval':True},
 'total_active_boxes':292849,'total_saved_bbox_visits':enum['total_saved_bbox_visits'],
 'all_candidate_classifications_independently_replayed':False,'raw_join_conditional_scopes_preserved':True,
 'report_sha256':sha(report_path),'evidence_sha256':{str(p):sha(p) for p in [sweep_path,enum_path,enum_review_path,report_path,secant_path,prior_review_path]},
 'script_sha256':sha(Path(__file__)),'actual_command':[sys.executable,str(Path(__file__).resolve()),*map(str,[sweep_path,enum_path,enum_review_path,report_path,secant_path,prior_review_path,out])]}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'script_sha256':result['script_sha256'],'output_sha256':sha(out)}))
