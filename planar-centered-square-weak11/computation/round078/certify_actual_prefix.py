#!/usr/bin/env python3
"""Single actual S20 input: complete Gamma prefix menu, directed exact areas."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,os,platform,sys,time
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ROOT=Path(__file__).resolve().parents[2]
CORE=ROOT/'computation/round075/certify_three_atom_collapse.py'
assert sha(CORE)=='a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456'
spec=importlib.util.spec_from_file_location('frozen_r075',CORE);area=importlib.util.module_from_spec(spec);spec.loader.exec_module(area)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);ap.add_argument('--bits',type=int,default=72);a=ap.parse_args();t=time.monotonic();ip=Path(a.input);data=json.loads(ip.read_text());points=[tuple(map(F,p)) for p in data['positions']];masses=list(map(F,data['masses']))
 assert F(data['threshold'])==1 and points[:400]==[(F(21*i,400),F(21*j,400)) for i in range(20) for j in range(20)] and masses[:400]==[F(1,400)]*400
 assert points[400:]==[(F(-71,200),F(0)),(F(0),F(-171,200))] and masses[400:]==[F(1,4),F(1,2)] and data['S_indices']==list(range(400))
 L=400*(1<<a.bits);unit=L//400;gamma=[-200*unit,-unit,-200*unit,-unit];rows=[]
 def clip(r):
  z=[max(r[0],gamma[0]),min(r[1],gamma[1]),max(r[2],gamma[2]),min(r[3],gamma[3])]
  return z if z[0]<z[1] and z[2]<z[3] else None
 def row(k,l,mask):
  old=k*l;coords=[]
  if old:coords=[(0,0),(21*(k-1),21*(l-1))]
  if mask&1:coords.append((-142,0))
  if mask&2:coords.append((0,-342))
  xmin=min(x for x,y in coords);xmax=max(x for x,y in coords);ymin=min(y for x,y in coords);ymax=max(y for x,y in coords)
  mass=F(old,400)+F(bool(mask&1),4)+F(bool(mask&2),2);span=F(max(xmax-xmin,ymax-ymin),400)
  td=max(xmax-xmin,ymax-ymin,2*(xmax+1),2*(ymax+1),2*(-200-xmin),2*(-200-ymin),0)
  active=mass>F(td*td,160000);r={'id':f'old{k}x{l}-outside{mask}','k':k,'l':l,'outside_mask':mask,'old_count':old,'mass':str(mass),'bbox_units_400':[xmin,xmax,ymin,ymax],'eligible':mass>span*span,'Gamma_strict_active':active,'twice_min_distance_units_400':td}
  assert mass<=1
  if r['eligible']:
   rad=mass*L*L;assert rad.denominator==1;n=math.isqrt(rad.numerator//4);up=n if 4*n*n==rad else n+1;assert 4*n*n<=rad<=4*up*up and up<=L//2
   r.update(radius_lower_units=n,radius_upper_units=up,inner=clip([xmax*unit-n,xmin*unit+n,ymax*unit-n,ymin*unit+n]),outer=clip([xmax*unit-up,xmin*unit+up,ymax*unit-up,ymin*unit+up]))
   if active:assert r['inner'] is not None
  else:assert not active
  return r
 for k in range(1,11):
  for l in range(1,11):
   for mask in range(4):rows.append(row(k,l,mask))
 residual=[row(0,0,mask) for mask in (1,2,3)];rows+=residual
 def rr(source,key):return [r[key] for r in source if r.get(key) is not None]
 fi,fo,ri,ro=rr(rows,'inner'),rr(rows,'outer'),rr(residual,'inner'),rr(residual,'outer')
 def enc(fullinner,fullouter):
  lo=area.event_area(fullinner,ro);hi=area.event_area(fullouter,ri)
  lo2,clo=area.midpoint_area(fullinner,ro);hi2,chi=area.midpoint_area(fullouter,ri)
  assert (lo,hi)==(lo2,hi2) and 0<=lo<=hi
  return {'lower':str(F(lo,L*L)),'upper':str(F(hi,L*L)),'lower_decimal':float(F(lo,L*L)),'upper_decimal':float(F(hi,L*L)),'width':str(F(hi-lo,L*L)),'two_algorithms_identical':True,'midpoint_cells':[clo,chi]}
 total=enc(fi,fo);joint=[]
 for r in rows:
  if r['old_count'] and r['outside_mask']==3 and r['Gamma_strict_active']:
   z=enc(rr([r],'inner'),rr([r],'outer'));joint.append({'id':r['id'],'k':r['k'],'l':r['l'],'area_outside_complete_residual':z,'classification':'positive_surviving_area' if F(z['lower'])>0 else 'fully_residual_covered_up_to_null_boundary' if F(z['upper'])==0 else 'unresolved_at_directed_precision'})
 pureold=[r for r in rows if r['old_count'] and not r['outside_mask']];oneoutside=[r for r in rows if r['old_count'] and r['outside_mask'] in (1,2)];alljoint=[r for r in rows if r['old_count'] and r['outside_mask']==3]
 groups={name:enc(rr(src,'inner'),rr(src,'outer')) for name,src in [('pure_old',pureold),('old_plus_one_outside',oneoutside),('old_plus_both_outside',alljoint)]}
 target=F(1,8);main=F(80003,640000);out={'schema':'r078-actual-S20-Gamma-prefix-v1','pid':os.getpid(),'python':sys.version,'platform':platform.platform(),'input_sha256':sha(ip),'source_sha256':sha(Path(__file__)),'frozen_area_core_sha256':sha(CORE),'bits':a.bits,'integer_scale':L,'Gamma':['-1/2','-1/400','-1/2','-1/400'],'full_menu_count':len(rows),'old_prefix_labels':400,'pure_outside_labels':3,'residual_label_count':3,'menu':rows,'full_inner_rectangles':fi,'full_outer_rectangles':fo,'residual_inner_rectangles':ri,'residual_outer_rectangles':ro,'active_counts_by_outside_mask':{str(mask):sum(r['Gamma_strict_active'] for r in rows[:400] if r['outside_mask']==mask) for mask in range(4)},'actual_Gamma_deletion_area':total,'groups_individually_after_complete_residual':groups,'joint_active_label_count':len(joint),'joint_label_classification':joint,'half_internal_P_mass':str(target),'main_corner_threshold':str(main),'strict_below_half_P':F(total['upper'])<target,'half_P_minus_upper':str(target-F(total['upper'])),'strict_below_main_threshold':F(total['upper'])<main,'main_threshold_minus_upper':str(main-F(total['upper'])),'algorithm_A':'vertical x-event direct interval union/difference','algorithm_B':'Cartesian compressed midpoint cell coverage','scope':'Actual unchanged S20+two outside, this Gamma only; complete menu uses r<1/2 justified by prior all-square caps. No collapsed old mass. No parameter search or general C2 improvement. Menu completeness proof is separately owned.','seconds':time.monotonic()-t}
 Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['pid','full_menu_count','active_counts_by_outside_mask','actual_Gamma_deletion_area','groups_individually_after_complete_residual','joint_active_label_count','joint_label_classification','strict_below_half_P','strict_below_main_threshold','seconds']},indent=2))
if __name__=='__main__':main()
