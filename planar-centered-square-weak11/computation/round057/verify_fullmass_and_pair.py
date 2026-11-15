#!/usr/bin/env python3
"""R57 source-pin whole-mass candidate premise and summarize exact derivative jumps."""
import argparse,gzip,hashlib,json,os,sys,platform
from fractions import Fraction as F
from decimal import Decimal,localcontext
from pathlib import Path
from full_area_sided_sweep import expr,serial

def load(path):
 raw=path.read_bytes();payload=gzip.decompress(raw) if path.suffix=='.gz' else raw
 return json.loads(payload),{'path':str(path.resolve()),'file_sha256':hashlib.sha256(raw).hexdigest(),'payload_sha256':hashlib.sha256(payload).hexdigest(),'bytes':len(raw)}
def main():
 a=argparse.ArgumentParser();a.add_argument('--joint',required=True);a.add_argument('--geometry',required=True);a.add_argument('--output',required=True);v=a.parse_args();here=Path(__file__).resolve().parent;g=Path(v.geometry)
 joint,jpin=load(Path(v.joint));ds=[];pins=[jpin]
 for name in ['one_axis_t0','one_axis_t1']:
  p=g/name/'slab.json.gz';d,pin=load(p);ds.append(d);pins.append(pin)
 a0,a1=ds
 for key in ['x_classes','y_classes','domain','s0','periodic','source_order_sha256']:assert a0[key]==a1[key]
 assert all({k:v for k,v in q.items() if k!='w'}=={k:v for k,v in r.items() if k!='w'} for q,r in zip(a0['source_catalog'],a1['source_catalog'],strict=True))
 maps=[{tuple(r['bbox']):box for r,box in zip(d['box_records'],d['boxes'],strict=True)} for d in ds]
 assert maps[0].keys()==maps[1].keys()=={tuple(r['bbox']) for r in joint['box_records']}
 assert all(F(d['H'])>0 and F(d['delta_empty'])>0 and d['status']=='complete' for d in ds)
 assert all(F(d['tail']['radius'])==F(joint['R']) and F(d['tail']['bound'])<F(d['tail']['lambda'])==F(joint['lambda']) for d in ds)
 t0=F(joint['t0']);R2=F(joint['R'])**2;leastq=None;leastgap=None
 for r,box in zip(joint['box_records'],joint['boxes'],strict=True):
  key=tuple(r['bbox']);q,b=F(r['q0']),F(r['q1']);qs=[q-b*t0,q+b*(1-t0)]
  assert 0<min(qs)<=max(qs)<R2
  leastq=min(qs) if leastq is None else min(leastq,*qs);leastgap=R2-max(qs) if leastgap is None else min(leastgap,R2-max(qs))
  for t,source in zip([F(0),F(1)],maps):
   expected=[expr([e[0],e[1],e[2],F(e[3])+F(e[4])*(t-t0),0]) for e in box]
   actual=[expr([*e,0]) for e in source[key]]
   assert expected==actual
 # All boxes omitted at both endpoints remain omitted for every t: their
 # squared activity threshold at fixed z is fixed, whereas bbox mass is affine.
 # A min with the fixed R cap does not change the implication. This is the
 # endpoint-product lemma from the pinned R52 whole-product proof below.
 proof=g/'wall_local_2d.md';audit=g/'audit_wall_2d.py';pins += [{'path':str(p.resolve()),'file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in [proof,audit]]
 M=F(joint['uniform_mass']['center_total_mass']);assert F(joint['uniform_mass']['total_derivative'])==0;lam=F(joint['lambda'])
 results={};sweep_pins=[]
 for side in ['left','right']:
  out,pin=load(here/f'full-area-{side}-256.json');results[side]=out;sweep_pins.append(pin)
  old,_=load(here/f'full-area-{side}.json')
  for key in ['A_wall_bounds','A_t_wall_bounds','beta_wall_bounds']:
   assert F(old[key][0])<=F(out[key][0])<=F(out[key][1])<=F(old[key][1])
  assert out['status']=='passed_certified_one_sided_limit'
 def jump(key):
  l,r=results['left'][key],results['right'][key];return F(r[0])-F(l[1]),F(r[1])-F(l[0])
 jumps={key:jump(key) for key in ['A_wall_bounds','A_t_wall_bounds','beta_wall_bounds']}
 assert jumps['A_wall_bounds'][0]<=0<=jumps['A_wall_bounds'][1]
 assert jumps['A_t_wall_bounds'][1]<0 and jumps['beta_wall_bounds'][0]>0
 assert all(F(r['A_t_wall_bounds'][0])>0 for r in results.values())
 out={'status':'passed','scope':'full union one-sided derivative kink for the specified complete periodic family; no new weak-type bound or extremizer','whole_mass_candidate_premise':{'endpoint_product_lemma':'fixed geometry plus affine bbox mass and identical endpoint active keys excludes every omitted bbox for all t in [0,1]; endpoint source/empty position radii intersect positively','active_bbox_keys':len(maps[0]),'all_candidate_endpoint_expressions_match_joint_affine_mass_model':True,'all_candidates_strictly_positive_root_branch_on_t_0_1':True,'minimum_q_on_t_0_1':leastq,'minimum_R_squared_minus_q_on_t_0_1':leastgap,'common_z0':a0['s0'],'common_position_radius':min(F(d[k]) for d in ds for k in ['H','delta_empty']),'uniform_fullscale_tail_bound_endpoint_pins':True,'source_pins':pins},'full_area_sweep_pins':sweep_pins,'A_t_left_bounds':results['left']['A_t_wall_bounds'],'A_t_right_bounds':results['right']['A_t_wall_bounds'],'right_minus_left_A_t_bounds':jumps['A_t_wall_bounds'],'right_minus_left_beta_bounds':jumps['beta_wall_bounds'],'area_limit_difference_bounds':jumps['A_wall_bounds'],'mass':M,'lambda':lam,'mass_t_derivative':'0','right_minus_left_J_t_bounds':[v*lam/M for v in jumps['A_t_wall_bounds']],'strict_conclusions':['A_t right is strictly smaller than A_t left','beta=A_z right is strictly larger than beta left','both one-sided A_t values remain positive','192-bit provisional enclosures contain final 256-bit enclosures'],'wall_equation':results['left']['wall_equation'],'wall_t_bounds':results['left']['wall_t_bounds'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'platform':platform.platform(),'actual_command':[sys.executable,*sys.argv]}
 with localcontext() as c:
  c.prec=35;dec=lambda q:str(Decimal(F(q).numerator)/Decimal(F(q).denominator))
  out['decimal_display_only']={key:[dec(q) for q in out[key]] for key in ['A_t_left_bounds','A_t_right_bounds','right_minus_left_A_t_bounds','right_minus_left_beta_bounds','right_minus_left_J_t_bounds','wall_t_bounds']}
 assert os.path.ismount('/Volumes/data');Path(v.output).write_text(json.dumps(serial(out),indent=2)+'\n');print(json.dumps(serial({key:out[key] for key in ['status','strict_conclusions','decimal_display_only','script_sha256']}),indent=2))
if __name__=='__main__':main()
