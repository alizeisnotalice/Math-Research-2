#!/usr/bin/env python3
"""One <=1000point n3 screen: exact squarecaps, floating actual Gamma area.

Local3-source full and residual catalogs retain7 and3 nonempty subsets.
All globalS3 bbox caps are checked exactly for each midpoint. NoB&B here.
"""
import argparse,hashlib,json,math,os,platform,sys,time
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def masks(points):
 axes=[]
 for axis in (0,1):
  xx=sorted(set(p[axis] for p in points));aa=[]
  for i,a in enumerate(xx):
   for b in xx[i:]:aa.append(sum(1<<j for j,p in enumerate(points) if a<=p[axis]<=b))
  axes.append(aa)
 return sorted({a&b for a in axes[0] for b in axes[1]}-{0})
def merge(iv):
 out=[]
 for a,b in sorted(iv):
  if not out or a>out[-1][1]:out.append([a,b])
  elif b>out[-1][1]:out[-1][1]=b
 return out
def intersect(a,b):
 i=j=0;z=0.
 while i<len(a) and j<len(b):
  z+=max(0.,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
  if a[i][1]<=b[j][1]:i+=1
  else:j+=1
 return z
def difference_area(full,res):
 events=defaultdict(list)
 for ch,rr in enumerate((full,res)):
  for i,(a,b,c,d) in enumerate(rr):events[a].append((ch,i,(c,d)));events[b].append((ch,i,None))
 if not events:return 0.
 active=[{},{}];last=min(events);area=0.
 for x,ee in sorted(events.items()):
  a=merge(active[0].values());b=merge(active[1].values());area+=(x-last)*(sum(t-s for s,t in a)-intersect(a,b))
  for ch,i,iv in ee:
   if iv is None:del active[ch][i]
   else:active[ch][i]=iv
  last=x
 assert not active[0] and not active[1];return area
def B(a):return min(max(F(1),(a+(k-1)*F(4,9))**2)-F(k*k,9) for k in (1,2,3))/2
def area(a,b):
 P=[(F(0),F(0)),(-a,F(0)),(F(0),-a)];M=[F(1,9),b,b];full=[];res=[]
 for mask in range(1,8):
  ii=[i for i in range(3) if mask>>i&1];mass=sum(M[i] for i in ii);xs=[P[i][0] for i in ii];ys=[P[i][1] for i in ii];span=max(max(xs)-min(xs),max(ys)-min(ys))
  if mass<=span*span:continue
  r=math.sqrt(float(mass))/2;rect=(max(-.5,float(max(xs))-r),min(-1/9,float(min(xs))+r),max(-.5,float(max(ys))-r),min(-1/9,float(min(ys))+r))
  if rect[0]<rect[1] and rect[2]<rect[3]:
   full.append(rect)
   if not mask&1:res.append(rect)
 return difference_area(full,res)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();start=time.monotonic();ip=Path(a.input);spec=json.loads(ip.read_text());assert spec['n']==3;Na=spec['a_midpoints'];Nt=spec['relative_b_midpoints'];assert Na*Nt<=1000;alo,ahi=map(F,spec['a_domain']);threshold=F(spec['comparison_threshold']);seed_points=[(F(4*i,9),F(4*j,9)) for i in range(3) for j in range(3)]+[(F(-1,2),F(0)),(F(0),F(-1,2))];closure_masks=masks(seed_points);rows=[]
 for i in range(Na):
  aa=alo+(ahi-alo)*F(2*i+1,2*Na);cap=B(aa);assert cap>0;P=[(F(4*k,9),F(4*j,9)) for k in range(3) for j in range(3)]+[(-aa,F(0)),(F(0),-aa)];boxdata=[]
  assert masks(P)==closure_masks
  for mask in closure_masks:
   ii=[k for k in range(11) if mask>>k&1];span=max(max(P[k][axis] for k in ii)-min(P[k][axis] for k in ii) for axis in (0,1));boxdata.append((mask,F(sum(k<9 for k in ii),9),sum(k>=9 for k in ii),max(F(1),span*span)))
  for j in range(Nt):
   t=F(2*j+1,2*Nt);bb=cap*t;slacks=[bound-internal-number*bb for mask,internal,number,bound in boxdata];assert min(slacks)>=0;value=area(aa,bb);rows.append({'a_index':i,'relative_b_index':j,'a':str(aa),'b':str(bb),'relative_b':str(t),'B3_a':str(cap),'area_float':value,'margin_to_1over18_float':float(threshold)-value,'all_global_caps_exact_checked':len(boxdata),'minimum_global_cap_slack':str(min(slacks))})
 top=sorted(rows,key=lambda r:r['area_float'],reverse=True)[:10];baseline=F(1,324);assert abs(area(F(1,9),F(0))-float(baseline))<1e-14
 out={'schema':'r075-orthogonal-n3-midpoint-screen-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'input_path':str(ip),'input_sha256':sha(ip),'parameters':spec,'count':len(rows),'complete_global_bbox_cap_masks':len(closure_masks),'total_exact_cap_checks':len(rows)*len(closure_masks),'full_local_subsets_per_point':7,'residual_local_subsets_per_point':3,'zero_outside_baseline_exact':'1/324','top10':top,'all_results':rows,'maximum_area_float':top[0]['area_float'],'minimum_margin_float':top[0]['margin_to_1over18_float'],'count_float_above_threshold':sum(r['area_float']>float(threshold) for r in rows),'seconds':time.monotonic()-start,'scope':'One finite961-midpoint screen. All shared fullS3 squarecaps exact; local complete full/residual Gamma deletion areas floating only. b is sampled throughout its cap interval, not assumed monotone. No certified maximum, full parameter proof, new bound or B&B admission claimed. a<=1/9 has B=0, while a=1 and b endpoints are not numerically sampled.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'pid':os.getpid(),'count':len(rows),'max_area_float':out['maximum_area_float'],'min_margin_float':out['minimum_margin_float'],'top':top[0],'seconds':out['seconds']},indent=2))
if __name__=='__main__':main()
