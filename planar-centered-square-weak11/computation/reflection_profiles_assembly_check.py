#!/usr/bin/env python3
"""Independent exact interval assembly of the completed R35/R41 edge/corner profiles."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys
assert os.path.ismount('/Volumes/data')
here=Path(__file__).resolve();ap=argparse.ArgumentParser();ap.add_argument('--tag',choices=['r035','r041'],default='r041');ap.add_argument('--bulk-source',type=Path);ap.add_argument('--profile-dir',type=Path,required=True);ap.add_argument('--axis-error',type=Path,required=True);ap.add_argument('--output',type=Path);a=ap.parse_args()
if a.output is None:a.output=here.parent/('results/reflection-edge-kernel-independent/'+a.tag+'-assembly.json')
input_hash={'r035':'5021c54704af45b3ef8289c67a245926473434640372a934b2f15c4485f0c026','r041':'ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3'}[a.tag]
edge_hash='ea664c7f9520b9983b13383fd032ff3110a6379b7dea381108ed49c88db9fede';corner_hash='53fe05c311c6601d4a7221e7518e3185b0b966ed737cd2a11c9916e5ce375fc3';union_hash='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'
states=['bulk','left','right','bottom','top','seam-old','seam-new','corner-lb','corner-lt','corner-rb','corner-rt','seam-bottom-old','seam-bottom-new','seam-top-old','seam-top-new'];records=[];z={};base=None
sp=a.profile_dir/(a.tag+'-seam.json');sb=sp.read_bytes();seam=json.loads(sb);seam_hash=hashlib.sha256(sb).hexdigest();assert seam['status']=='PASS' and seam['input_sha256']==input_hash
for state in states:
    reused=(state=='bulk' and a.bulk_source is not None)
    path=a.bulk_source if reused else a.profile_dir/(a.tag+'-'+state+'.json');blob=path.read_bytes();o=json.loads(blob)
    if reused:
        assert o['input_sha256']==input_hash and o['code_sha256']==union_hash and F(o['level'])==1
        params=(F(o['period_x']),F(o['period_y']),F(o['mass']))
        win=(F(0),params[0],F(0),params[1])
    else:
        expected=corner_hash if state.startswith('corner-') or state.startswith('seam-bottom-') or state.startswith('seam-top-') else edge_hash
        assert o['status']=='PASS' and o['state']==state and o['input_sha256']==input_hash
        assert o['script_sha256']==expected and o['union_module_sha256']==union_hash and F(o['lam'])==1
        if state.startswith('seam-'):assert o['seam_sha256']==seam_hash and F(o['rational_axis'])==F(seam['rational_reflection_axis'])
        params=tuple(F(o[k]) for k in ['Px','Py','cell_mass']);win=list(map(F,o['window']))
    lo,hi=F(o['inner_area']),F(o['outer_area']);assert 0<=lo<=hi<=(win[1]-win[0])*(win[3]-win[2])
    assert F(o['tail_radius'])==(15 if state.endswith('-new') else 10) and F(o['tail_bound'])<1
    if base is None:base=params
    else:assert base==params
    z[state]=(lo,hi);records.append(dict(path=str(path.resolve()),bytes=len(blob),sha256=hashlib.sha256(blob).hexdigest(),state=state,reused_frozen_bulk=reused))
A=z['bulk'];E=(sum(z[s][0] for s in ['left','right','bottom','top']),sum(z[s][1] for s in ['left','right','bottom','top']))
Cc=(sum(z[s][0] for s in ['corner-lb','corner-lt','corner-rb','corner-rt']),sum(z[s][1] for s in ['corner-lb','corner-lt','corner-rb','corner-rt']))
b=(E[0]-8*A[1],E[1]-8*A[0]);g=(z['seam-new'][0]-z['seam-old'][1],z['seam-new'][1]-z['seam-old'][0])
G=(sum(z['seam-'+t+'-new'][0]-z['seam-'+t+'-old'][1] for t in ['bottom','top']),sum(z['seam-'+t+'-new'][1]-z['seam-'+t+'-old'][0] for t in ['bottom','top']))
ab=a.axis_error.read_bytes();axis=json.loads(ab);assert axis['status']=='PASS' and axis['input_sha256']==input_hash and axis['seam_sha256']==seam_hash and axis['edge_state_code_sha256']==edge_hash and axis['corner_state_code_sha256']==corner_hash
assert F(axis['delta_upper'])<=F(1,1<<64)
eg,eG=F(axis['gamma_change_upper']),F(axis['Gamma_change_upper'])
# Recheck every complete-catalog count and per-bbox endpoint bound from the saved exact rows.
by={}
for row in axis['states']:
    nx,ny=row['nx'],row['ny'];K=nx*(nx+1)//2*ny*(ny+1)//2;assert K==row['all_bbox_upper_count']
    bound=2*F(axis['delta_upper'])*F(row['height_bound'])*K;assert bound==F(row['area_change_upper']);by[row['state']]=bound
assert eg==by['seam-old']+by['seam-new'] and eG==sum(v for k,v in by.items() if k not in ('seam-old','seam-new'))
slope=(b[0]+g[0],b[1]+g[1]);slope_true=(slope[0]-eg,slope[1]+eg)
# At N=6 all gamma terms cancel exactly: F6=2E-32A+Cc+Gamma.
f6=(2*E[0]-32*A[1]+Cc[0]+G[0],2*E[1]-32*A[0]+Cc[1]+G[1]);f6_true=(f6[0]-eG,f6[1]+eG)
assert slope_true[1]<-9 and f6_true[1]<-53
px,py,mass=base;assert 2*px>10 and 2*py>10 and 3*py>15 and 5*px+F(seam['c_edge_lower'])>30
rec=lambda v:dict(lower=str(v[0]),upper=str(v[1]),lower_float=float(v[0]),upper_float=float(v[1]))
result=dict(status='PASS',python=sys.version,arithmetic='Fraction exact interval arithmetic',random_seed=None,floating_tolerance=None,input_sha256=input_hash,tag=a.tag,source_profiles=records,seam_source=dict(path=str(sp.resolve()),bytes=len(sb),sha256=seam_hash),axis_error_source=dict(path=str(a.axis_error.resolve()),bytes=len(ab),sha256=hashlib.sha256(ab).hexdigest()),A=rec(A),E=rec(E),Ccorner=rec(Cc),b=rec(b),gamma_rational_axis=rec(g),Gamma_rational_axis=rec(G),slope_rational_axis=rec(slope),slope_true_axis=rec(slope_true),F6_rational_axis=rec(f6),F6_true_axis=rec(f6_true),gamma_axis_error=str(eg),Gamma_axis_error=str(eG),clean_rational_bounds={'leading_upper_strict':'-9','F6_upper_strict':'-53','ratio_deficit_strict_lower':'(9*N-1)/(mass_per_cell*N^2)'},conclusion='For every integer N>=6, F_N(true C)<-(9*N-1). The reflected finite-block ratio is below its exact old periodic bulk ratio by more than (9*N-1)/(mass_per_cell*N^2).',scope='Only the specified fixed input, literal one-step right reflection at its true or certified outward axis, and N>=6. N=1..5, other templates, other reflection strategies, and general C2 are not decided.',script_sha256=hashlib.sha256(here.read_bytes()).hexdigest())
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','slope_true_axis','F6_true_axis','conclusion']},ensure_ascii=False,indent=2))
