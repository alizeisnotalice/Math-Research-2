#!/usr/bin/env python3
"""Exact full-segment window geometry and constants, not a Borel area proof.

Captured-index lists refer only to the four old points and eight named
segment families. Additional rho may enter capacity windows and only
strengthens their necessary mass constraints. Its exclusion from corner
witnesses is an analytic support-window statement in the external proof.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys,time
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class A:
 def __init__(self,const=0,coeff=None):self.const=F(const);self.coeff={k:F(v) for k,v in (coeff or {}).items() if v}
 def __add__(self,other):
  if not isinstance(other,A):other=A(other)
  d=self.coeff.copy()
  for k,v in other.coeff.items():d[k]=d.get(k,F(0))+v
  return A(self.const+other.const,d)
 __radd__=__add__
 def __neg__(self):return A(-self.const,{k:-v for k,v in self.coeff.items()})
 def __sub__(self,other):return self+-other if isinstance(other,A) else self+A(-F(other))
 def __rsub__(self,other):return A(other)+-self
 def interval(self,ranges):
  lo=hi=self.const
  for k,v in self.coeff.items():
   a,b=ranges[k];lo+=v*(a if v>0 else b);hi+=v*(b if v>0 else a)
  return lo,hi
 def serial(self):return {'constant':str(self.const),'coefficients':{k:str(v) for k,v in sorted(self.coeff.items())}}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--formal-proof',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();started=time.monotonic();ip=Path(a.input);data=json.loads(ip.read_text());assert data['schema']=='r082-robust-shared-axis-geometry-v1';proof=Path(a.formal_proof);assert proof.is_file()
 ranges={'w':tuple(map(F,data['w_closed_verification_range']))}
 for i in range(4):ranges[f'd{i}']=tuple(map(F,data['near_distance_range']));ranges[f'D{i}']=tuple(map(F,data['far_distance_range']))
 w=A(coeff={'w':1});ds=[A(coeff={f'd{i}':1}) for i in range(4)];Ds=[A(coeff={f'D{i}':1}) for i in range(4)];zero=A(0)
 def rotate_point(p):x,y=p;return w-y,x
 def rotate_box(b):a,b,c,d=b;return w-d,w-c,a,b
 def rotated(obj,count,fn):
  for _ in range(count):obj=fn(obj)
  return obj
 points=[(zero,zero),(w,zero),(w,w),(zero,w)]
 for i in range(4):points.extend([rotated((-ds[i],zero),i,rotate_point),rotated((zero,-Ds[i]),i,rotate_point)])
 def audit(box,expected,kind,label):
  aa,bb,cc,dd=box;records=[]
  for j,(x,y) in enumerate(points):
   inequalities=[x-aa,bb-x,y-cc,dd-y];bounds=[f.interval(ranges) for f in inequalities]
   if j in expected:
    assert all(lo>0 if kind=='open' else lo>=0 for lo,hi in bounds);reason='uniform_inside'
   else:
    facets=[k for k,(lo,hi) in enumerate(bounds) if hi<0];assert facets;reason=f'uniform_outside_facet_{facets[0]}'
   records.append({'source_index':j,'expected_captured':j in expected,'reason':reason,'four_signed_inclusion_facets':[{'affine':f.serial(),'minimum':str(lo),'maximum':str(hi)} for f,(lo,hi) in zip(inequalities,bounds)]})
  return {'label':label,'boundary':kind,'box_affine':[z.serial() for z in box],'expected_capture_indices':sorted(expected),'source_checks':records,'all_affine_extrema_verified':True}
 windows=[];weight_counts=[0]*8;dmax=ranges['d0'][1];Dmax=ranges['D0'][1]
 for i in range(4):
  # A single fixed window must contain the WHOLE support segments.
  # It may depend on the old width w, never on a sampled source d_i/D_i.
  prev=(i-1)%4;near=4+2*i;farprev=5+2*prev;side=rotated((A(-Dmax),A(1-Dmax),zero,A(1)),i,rotate_box);expected={i,prev,near,farprev};assert all(set(z.coeff)<=set(['w']) for z in side);windows.append(audit(side,expected,'closed',f'side_{i}'));weight_counts[2*i]+=1;weight_counts[2*prev+1]+=1
  local=rotated((A(-dmax),A(1-dmax),A(-Dmax),A(1-Dmax)),i,rotate_box);assert all(set(z.coeff)<=set(['w']) for z in local);windows.append(audit(local,{i,near,near+1},'closed',f'corner_unit_{i}'))
  P=rotated((A(-1),w-F(1,2),A(-1),w-F(1,2)),i,rotate_box);windows.append(audit(P,{i,near,near+1},'open',f'effective_source_window_{i}'))
 assert weight_counts==[1]*8
 m=F(data['old_corner_mass']);amax=F(data['near_mass_upper_from_side_caps']);assert m==F(1,4) and amax==F(1,2)
 # u^2=(m+a)/4<=3/16 and z^2=a/4<=1/8. Positive rational
 # bounds u<7/16,z<3/8 imply u+2z<19/16<6/5.
 ub,zb=F(7,16),F(3,8);assert ub*ub>(m+amax)/4 and zb*zb>amax/4;lip=ub+2*zb;assert lip==F(19,16)<F(6,5)
 e0=F(data['base_crop_e']);crop=F(1,4)-(F(1,2)-e0)**2;assert crop==F(999,1000000)<F(1,1000)
 d0=F(data['base_near_distance']);displace=F(data['near_distance_displacement_bound']);assert max(abs(ranges['d0'][0]-d0),abs(ranges['d0'][1]-d0))==displace
 nearbound=F(data['base_near_area_upper_from_separate_proof'])+crop+F(6,5)*displace;assert nearbound==F(74499,1000000)<F(149,2000)<F(3,40)
 # Far-only actual deletion is bounded by its mixed full rectangle.
 sqrt3lower=F(433,250);assert sqrt3lower*sqrt3lower<3;farbound=F(3,8)-F(17,80)*sqrt3lower;assert farbound==F(139,20000)<F(7,1000)
 theta=ranges['D0'][0]**2-m;kappa=F(1,2)-ranges['D0'][0]/4;assert theta==F(189,400) and kappa==F(23,80)
 hinge=[max(F(0),min(F(2),F(3*k,4))-k*theta) for k in range(5)];assert hinge==[F(0),F(111,400),F(111,200),F(233,400),F(11,100)]
 joint=kappa*max(hinge);assert joint==F(5359,32000)<F(21,125)
 total=4*F(3,40)+4*F(7,1000)+joint;assert total==F(3171,6400)<F(62,125)<F(1,2)
 out={'schema':'r082-robust-geometry-constants-v1','pid':os.getpid(),'python':sys.version,'source_sha256':sha(__file__),'input_sha256':sha(ip),'formal_proof_path':str(proof),'formal_proof_sha256':sha(proof),'parameter_ranges':{k:list(map(str,v)) for k,v in ranges.items()},'actual_w_open_interval':'3/4<w<1; closed range used only for safe affine extrema','windows':windows,'window_count':12,'source_window_incidence_checks':144,'signed_affine_facet_checks':576,'side_weight_multiplicity':weight_counts,'necessary_weight_constraints':{'side':'a_i+b_(i-1)<=1/2','corner':'a_i+b_i<=3/4','summed_side':'sum_i(a_i+b_i)<=2','nonnegative_weights':True,'sufficient_for_global_caps':False},'robust_constants':{'u_rational_upper':str(ub),'z_rational_upper':str(zb),'distance_lipschitz_upper':str(lip),'coarser_lipschitz_used':'6/5','crop_area_loss':str(crop),'distance_displacement':str(displace),'base_near_bound_assumed_from_separate_proof':data['base_near_area_upper_from_separate_proof'],'near_bound_after_perturbation':str(nearbound),'far_full_rectangle_upper':str(farbound),'joint_threshold':str(theta),'joint_kappa_upper':str(kappa),'active_count_hinge_bounds':list(map(str,hinge)),'joint_sum_upper':str(joint),'sum_using_075_near_and_007_far':str(total),'coarser_aggregate_upper':'62/125','aggregate_strict_below_half':True},'all_checks_pass':True,'area_parameter_samples':0,'branch_and_bound_nodes':0,'random_seed':None,'numerical_tolerance':None,'seconds':time.monotonic()-started,'scope':'Independent affine containment/exclusion for all robust geometry parameters and rational constant checks conditional on separate analytic near/far/joint lemmas. Does not duplicate the base near-only three-branch area proof, enumerate mass vectors or claim global cap sufficiency from these selected windows.'}
 out['window_capture_sets_exact_only_among_named_support_families']=True
 out['entire_segments_verified_by_affine_extrema_not_point_samples']=True
 out['capacity_windows_independent_of_all_segment_point_variables']=True
 out['additional_rho_scope']='Not enumerated. May enter capacity windows, strengthening the necessary inequalities by nonnegativity. No rho mass enters corner witnesses by the separate analytic P-window containment theorem.'
 out['geometry_check_is_not_Borel_area_proof']=True
 out['pre_execution_prepared_source_sha256']='d77f409b8582cbfa1b4e39cd94f287583aa425d95673222ee1c8b66d192e7e9a'
 op=Path(a.out);assert not op.exists();op.parent.mkdir(parents=True,exist_ok=True);op.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['pid','window_count','source_window_incidence_checks','signed_affine_facet_checks','necessary_weight_constraints','robust_constants','all_checks_pass','seconds']},indent=2))
if __name__=='__main__':main()
