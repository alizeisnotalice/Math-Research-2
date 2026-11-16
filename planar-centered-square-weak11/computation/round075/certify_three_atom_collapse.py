#!/usr/bin/env python3
"""One directed72bit three-atom counterexample, two exact area algorithms.

All7 nonempty full-source subsets and all3 residual subsets are retained.
No maximal-block squarecaps are imposed; this is only an uncapped auxiliary
collapse test. A second area algorithm enumerates compressed2D midpoint
cells, independently of the first vertical event/interval-union sweep.
"""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import argparse,hashlib,json,math,os,platform,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def merge(iv):
 out=[]
 for a,b in sorted(iv):
  if not out or a>out[-1][1]:out.append([a,b])
  elif b>out[-1][1]:out[-1][1]=b
 return out
def intersection(a,b):
 i=j=0;out=0
 while i<len(a) and j<len(b):
  out+=max(0,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
  if a[i][1]<=b[j][1]:i+=1
  else:j+=1
 return out
def event_area(full,res):
 events=defaultdict(list)
 for ch,rr in enumerate((full,res)):
  for j,(a,b,c,d) in enumerate(rr):events[a].append((ch,j,(c,d)));events[b].append((ch,j,None))
 active=[{},{}];last=min(events);area=0
 for x,ee in sorted(events.items()):
  if x>last:
   f=merge(active[0].values());r=merge(active[1].values());area+=(x-last)*(sum(b-a for a,b in f)-intersection(f,r))
  for ch,j,value in ee:
   if value is None:del active[ch][j]
   else:active[ch][j]=value
  last=x
 assert not active[0] and not active[1];return area
def midpoint_area(full,res):
 xx=sorted({q for rr in full+res for q in rr[:2]});yy=sorted({q for rr in full+res for q in rr[2:]});area=0;cells=0
 for a,b in zip(xx,xx[1:]):
  for c,d in zip(yy,yy[1:]):
   x2=a+b;y2=c+d;cells+=1
   inside=lambda rr:any(2*u<x2<2*v and 2*s<y2<2*t for u,v,s,t in rr)
   if inside(full) and not inside(res):area+=(b-a)*(d-c)
 return area,cells
def SW(rr):return [(a,min(b,0),c,min(d,0)) for a,b,c,d in rr if a<min(b,0) and c<min(d,0)]
def enc(lo,hi,D):
 assert 0<=lo<=hi;return {'lower':str(F(lo,D*D)),'upper':str(F(hi,D*D)),'lower_decimal':float(F(lo,D*D)),'upper_decimal':float(F(hi,D*D)),'width':str(F(hi-lo,D*D))}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);ap.add_argument('--bits',type=int,default=72);a=ap.parse_args();start=time.monotonic();ip=Path(a.input);data=json.loads(ip.read_text());points=[tuple(map(F,p)) for p in data['positions']];masses=list(map(F,data['masses']));assert points==[(F(0),F(0)),(F(-2),F(0)),(F(0),F(-2))] and masses==[1,4,4];assert F(data['threshold'])==1;D=1<<a.bits;catalog=[]
 for mask in range(1,8):
  jj=[j for j in range(3) if mask>>j&1];mass=sum(masses[j] for j in jj);box=(min(points[j][0] for j in jj),max(points[j][0] for j in jj),min(points[j][1] for j in jj),max(points[j][1] for j in jj));span=max(box[1]-box[0],box[3]-box[2]);eligible=mass>span*span;row={'mask':mask,'indices':jj,'mass':str(mass),'bbox':list(map(str,box)),'span_squared':str(span*span),'eligible':eligible}
  if eligible:
   radnum=mass*D*D;assert radnum.denominator==1;k=math.isqrt(int(radnum)//4);assert 4*k*k<=radnum<4*(k+1)*(k+1);upper=k if 4*k*k==radnum else k+1;a0,b,c,d=(int(v*D) for v in box);row.update(radius_lower_units=k,radius_upper_units=upper,inner=[b-k,a0+k,d-k,c+k],outer=[b-upper,a0+upper,d-upper,c+upper])
  catalog.append(row)
 assert len(catalog)==7 and all(r['eligible'] for r in catalog)
 fi=SW([r['inner'] for r in catalog]);fo=SW([r['outer'] for r in catalog]);ri=SW([r['inner'] for r in catalog if not r['mask']&1]);ro=SW([r['outer'] for r in catalog if not r['mask']&1]);lo=event_area(fi,ro);hi=event_area(fo,ri);lo2,cellslo=midpoint_area(fi,ro);hi2,cellshi=midpoint_area(fo,ri);assert (lo,hi)==(lo2,hi2);area=enc(lo,hi,D);assert F(area['lower'])>F(1,2)
 # Independent exact hand formula: Delta=4sqrt2+sqrt5/2-25/4.
 vlow=F(707,500);ulow=F(559,500);assert vlow*vlow<2 and ulow*ulow<F(5,4);simple_lower=4*vlow+ulow-F(25,4);assert simple_lower==F(131,250)>F(1,2)
 vk=math.isqrt(2*D*D);uk=math.isqrt(5*D*D//4);formula_lower=4*F(vk,D)+F(uk,D)-F(25,4);formula_upper=4*F(vk+1,D)+F(uk+1,D)-F(25,4);assert max(F(area['lower']),formula_lower)<=min(F(area['upper']),formula_upper)
 out={'schema':'r075-three-atom-uncapped-collapse-counterexample-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'input_path':str(ip),'input_sha256':sha(ip),'bits':a.bits,'coordinate_scale':D,'complete_full_subset_count':7,'complete_residual_subset_count':3,'full_subset_catalog':catalog,'SW_full_inner_rectangles':fi,'SW_full_outer_rectangles':fo,'SW_residual_inner_rectangles':ri,'SW_residual_outer_rectangles':ro,'algorithmA':'verticalxevent/directyinterval union and subtraction','algorithmB':'compressed2D Cartesian midpoint cells, exact integer coverage predicates','algorithmA_Delta_SW':area,'algorithmB_Delta_SW':enc(lo2,hi2,D),'algorithmB_cells_lower':cellslo,'algorithmB_cells_upper':cellshi,'algorithms_integer_areas_identical':True,'internal_mass':'1','auxiliary_half_mass':'1/2','strict_violation':True,'strict_excess_lower':str(F(area['lower'])-F(1,2)),'exact_hand_formula':'4sqrt(2)+sqrt(5)/2-25/4','hand_formula_bracket':{'lower':str(formula_lower),'upper':str(formula_upper)},'simple_rational_strict_lower':'131/250','simple_rational_strict_excess':'3/125','rational_root_checks':{'2_minus_707over500_squared':str(2-vlow*vlow),'5over4_minus_559over500_squared':str(F(5,4)-ulow*ulow)},'squarecap_countercheck':{'closed_square':['-2','0','-2','0'],'captured_total_mass':'9','square_area':'4','violates_global_maxblock_normalization':True},'actual_global_maximum_eligible_mass':'9','S_is_global_maximum':False,'new_random_batch_inputs':0,'directed_family_inputs_evaluated':1,'seconds':time.monotonic()-start,'scope':'Strict counterexample only to the uncapped singleinternal plus twooutside collapse auxiliary bound. Not an admissible maximaluniformgrid input, not a counterexample to3.5, and not a C2 bound. All full/residual subsets including triple retained; no selected-radius proxy.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'pid':os.getpid(),'seconds':out['seconds'],'Delta_SW':area,'half_mass_violation':True,'simple_strict_lower':'131/250','random_batch_started':False,'further_d_stopped':True},indent=2))
if __name__=='__main__':main()
