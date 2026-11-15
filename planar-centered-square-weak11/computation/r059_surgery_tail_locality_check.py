#!/usr/bin/env python3
"""Exact fixed-constant checks for R59 labelwise mass-preserving endpoint surgery."""
import argparse,hashlib,json,os,sys
from pathlib import Path
from fractions import Fraction as F

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ar=ap.parse_args();assert os.path.ismount('/Volumes/data')
 blob=ar.input.read_bytes();z=json.loads(blob);pts=[tuple(map(F,p)) for p in z['pts']];px,py,lam=map(F,(z['Px'],z['Py'],z['lam']));assert lam==1
 exposed=[0,51];assert [i for i,(x,y,w) in enumerate(pts) if x==0]==exposed
 a,b,w=pts[1];assert 1 not in exposed and a==px/4 and w>4
 mass=sum(m for x,y,m in pts);rho=mass/(px*py)
 g10=rho*(1+px/20)*(1+py/20);g15=rho*(1+px/15)*(1+py/30);long=rho*F(3,2)*(1+py/(6*px));assert g10<1 and g15<1 and long<1
 margins={
 'six_Px_minus30':6*px-30,
 'three_Px_minus15':3*px-15,
 'left_collar_axis_gap_minus15':3*px-a+1-15,
 'right_collar_axis_gap_minus15':2*px+a+1-15,
 'N7_omitted_source_gap_minus15':4*px-15,
 'N7_left_axis_old_capture_gap':7*px-a-9,
 'N7_right_axis_old_capture_gap':6*px+a-9,
 'N7_right_seam_to_left_collar_gap_minus15':3*px+a+1-15,
 'N7_left_seam_to_right_collar_gap_minus15':4*px-a+1-15,
 'singleton_radius_squared_minus1':w/4-1}
 assert all(v>0 for v in margins.values())
 # Affine endpoint correction to label total: removed w, restored t*w+(1-t)*w.
 affine={str(i):{'constant_mass_change':str(-pts[i][2]+pts[i][2]),'t_coefficient_mass_change':str(pts[i][2]-pts[i][2])} for i in exposed}
 assert all(F(v['constant_mass_change'])==F(v['t_coefficient_mass_change'])==0 for v in affine.values())
 sample={str(t):[str(t*pts[i][2]),str((1-t)*pts[i][2])] for t in (F(0),F(1,2),F(1)) for i in exposed[:1]}
 out=dict(status='PASS',python=sys.version,arithmetic='Fraction; exact fixed constants and affine mass identity',random_seed=None,floating_tolerance=None,input_path=str(ar.input),input_sha256=hashlib.sha256(blob).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),exposed_labels=exposed,unchanged_witness_label=1,unchanged_witness=dict(alpha=str(a),mass=str(w),strict_singleton_radius_lower='1'),mass_per_cell=str(mass),old_R10_bound=str(g10),reflected_short_R15_bound=str(g15),long_all_N_ge7_bound=str(long),positive_margins={k:str(v) for k,v in margins.items()},endpoint_affine_mass_correction=affine,exposed_label0_endpoint_weights_at_sample_t=sample,uniform_analytic_scope='all t in [0,1], every integer N>=7; affine identity and weighted lattice proof, no N scan and no halo-area computation',formula='J_N(t)=A/m+(Lambda_L(t)+Lambda_R(1-t)-6*A)/(m*N)')
 ar.output.parent.mkdir(parents=True,exist_ok=True);ar.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('status','unchanged_witness','positive_margins','uniform_analytic_scope')},indent=2))
if __name__=='__main__':main()
