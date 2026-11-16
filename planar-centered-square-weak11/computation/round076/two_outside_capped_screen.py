#!/usr/bin/env python3
"""One finite n3 actual-corner screen: complete exact caps, floating areas.

Stop at the first floating signal above1/18+1e-5. No rigorous area decision
is made here. All7full and3residual local labels share one mass vector.
"""
from fractions import Fraction as F
from collections import defaultdict,Counter
from pathlib import Path
import argparse,gzip,hashlib,itertools,json,math,os,platform,random,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(p,q):return tuple(sorted((tuple(p),tuple(q))))
def pairs(params):
 N=params['negative_levels_units'];T=params['arm_other_levels_units'];L=[(x,y) for x in N for y in T];B=[(y,x) for x in N for y in T];C=list(itertools.product(N,N));rng=random.Random(params['pair_seed']);out=[];seen=set()
 def put(category,pp):
  for p,q in pp:
   z=canon(p,q);assert p!=q and z not in seen;seen.add(z);out.append((category,z))
 symmetric=[((x,y),(y,x)) for x in N for y in T];put('symmetric_cross',symmetric);remaining=[(p,q) for p in L for q in B if canon(p,q) not in seen];put('additional_cross',rng.sample(remaining,params['pair_counts']['additional_cross']));put('same_left',rng.sample(list(itertools.combinations(L,2)),params['pair_counts']['same_left']));put('same_bottom',rng.sample(list(itertools.combinations(B,2)),params['pair_counts']['same_bottom']));put('corner_plus_arm',rng.sample([(p,q) for p in C for q in L+B],params['pair_counts']['corner_plus_arm']));put('two_corners',rng.sample(list(itertools.combinations(C,2)),params['pair_counts']['two_corners']));assert len(out)==params['position_pair_limit']==1024
 return out
def ids(mask):
 out=[]
 while mask:b=mask&-mask;out.append(b.bit_length()-1);mask-=b
 return out
def caps(P,D):
 axes=[]
 for axis in (0,1):
  groups=defaultdict(int)
  for j,p in enumerate(P):groups[p[axis]]|=1<<j
  xx=sorted(groups);aa=[]
  for i in range(len(xx)):
   m=0
   for x in xx[i:]:m|=groups[x];aa.append(m)
  axes.append(aa)
 rows=[];bounds={1:None,2:None,3:None}
 for mask in sorted({a&b for a in axes[0] for b in axes[1]}-{0}):
  ii=ids(mask);span=max(max(P[j][axis] for j in ii)-min(P[j][axis] for j in ii) for axis in (0,1));tag=(mask>>9)&3;fixed=(mask&511).bit_count();rhs=max(D,span*span)-D//9*fixed;assert rhs>=0;rows.append((tag,F(rhs,D)))
  if tag and (bounds[tag] is None or rhs<bounds[tag]):bounds[tag]=rhs
 assert all(v is not None for v in bounds.values());return rows,tuple(F(bounds[k],D) for k in (1,2,3))
def mass_points(V1,V2,U):
 cc=[(F(0),F(0)),(V1,F(0)),(F(0),V2),(V1,V2),(F(0),U),(U,F(0)),(V1,U-V1),(U-V2,V2)];vertices=sorted({(x,y) for x,y in cc if 0<=x<=V1 and 0<=y<=V2 and x+y<=U});assert vertices
 center=(sum(p[0] for p in vertices)/len(vertices),sum(p[1] for p in vertices)/len(vertices));hi1=max(vertices,key=lambda p:(p[0],sum(p)));hi2=max(vertices,key=lambda p:(p[1],sum(p)));interior=[center,(center[0]/2,center[1]/2),((center[0]+hi1[0])/2,(center[1]+hi1[1])/2),((center[0]+hi2[0])/2,(center[1]+hi2[1])/2)];allpoints=[]
 for label,pp in [('vertex',vertices),('internal',interior)]:
  for p in pp:
   if p not in [v for k,v in allpoints]:allpoints.append((label,p))
 assert len(allpoints)<=9;return vertices,allpoints
def merge(iv):
 out=[]
 for a,b in sorted(iv):
  if not out or a>out[-1][1]:out.append([a,b])
  elif b>out[-1][1]:out[-1][1]=b
 return out
def intersection(a,b):
 i=j=0;z=0.
 while i<len(a) and j<len(b):
  z+=max(0.,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
  if a[i][1]<=b[j][1]:i+=1
  else:j+=1
 return z
def area(full,res):
 events=defaultdict(list)
 for ch,rr in enumerate((full,res)):
  for j,(a,b,c,d) in enumerate(rr):events[a].append((ch,j,(c,d)));events[b].append((ch,j,None))
 if not events:return (0.,0.,0.)
 active=[{},{}];last=min(events);out=[0.,0.,0.]
 for x,ee in sorted(events.items()):
  f=merge(active[0].values());r=merge(active[1].values());lf=sum(b-a for a,b in f);lr=sum(b-a for a,b in r);width=x-last;out[0]+=width*lf;out[1]+=width*lr;out[2]+=width*(lf-intersection(f,r))
  for ch,j,v in ee:
   if v is None:del active[ch][j]
   else:active[ch][j]=v
  last=x
 assert not active[0] and not active[1];return tuple(out)
def triple_delta(P):
 a=min(p[0] for p in P);b=max(p[0] for p in P);c=min(p[1] for p in P);d=max(p[1] for p in P)
 return max((b-a)/2,(d-c)/2,b+F(1,9),d+F(1,9),F(-1,2)-a,F(-1,2)-c,F(0))
def actual_area(q1,q2,b1,b2):
 P=[(F(0),F(0)),q1,q2];M=[F(1,9),b1,b2];full=[];res=[]
 for mask in range(1,8):
  jj=[j for j in range(3) if mask>>j&1];mass=sum(M[j] for j in jj);xs=[P[j][0] for j in jj];ys=[P[j][1] for j in jj];span=max(max(xs)-min(xs),max(ys)-min(ys))
  if mass<=span*span:continue
  r=math.sqrt(float(mass))/2;rr=(max(-.5,float(max(xs))-r),min(-1/9,float(min(xs))+r),max(-.5,float(max(ys))-r),min(-1/9,float(min(ys))+r))
  if rr[0]<rr[1] and rr[2]<rr[3]:
   full.append(rr)
   if not mask&1:res.append(rr)
 return area(full,res)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);ap.add_argument('--data',required=True);a=ap.parse_args();start=time.monotonic();ip=Path(a.input);params=json.loads(ip.read_text());assert params['n']==3;C=params['coordinate_denominator'];D=C*C;assert D%9==0;S=[(C*4*i//9,C*4*j//9) for i in range(3) for j in range(3)];pool=pairs(params);threshold=F(params['auxiliary_threshold']);trigger=threshold+F(params['floating_trigger_excess']);cache={};records=[];top=[];counts=Counter();capchecks=0;triggered=None;masspoints_total=0
 for index,(category,(Q1,Q2)) in enumerate(pool):
  for p in (Q1,Q2):assert -900<p[0]<350 and -900<p[1]<350 and (p[0]<-100 or p[1]<-100) and p not in S
  assert Q1!=Q2;P=S+[Q1,Q2];rows,(V1,V2,U)=caps(P,D);vertices,mm=mass_points(V1,V2,U);q1=tuple(F(x,C) for x in Q1);q2=tuple(F(x,C) for x in Q2);delta=triple_delta([(F(0),F(0)),q1,q2]);rec={'pair_index':index,'category':category,'positions_units':[Q1,Q2],'V1':str(V1),'V2':str(V2),'U':str(U),'full_bbox_cap_count':len(rows),'mass_vertices':[[str(x),str(y)] for x,y in vertices],'triple_min_farthest_distance':str(delta),'all_feasible_mass_points_have_null_formal_triple':F(1,9)+min(U,V1+V2)<=4*delta*delta,'mass_points':[]}
  for kind,(b1,b2) in mm:
   labelmass={0:F(0),1:b1,2:b2,3:b1+b2};slacks=[rhs-labelmass[tag] for tag,rhs in rows];assert min(slacks)>=0;capchecks+=len(rows);masspoints_total+=1;active=F(1,9)+b1+b2>4*delta*delta
   if active:assert F(1,9)+b1+b2<=1
   key=tuple(sorted((p[0],p[1],str(b)) for p,b in [(Q1,b1),(Q2,b2)] if b>0));cached=key in cache
   if not cached:cache[key]=actual_area(q1,q2,b1,b2)
   full,res,value=cache[key];r={'kind':kind,'masses':[str(b1),str(b2)],'formal_triple_Gamma_active':active,'cached_actual_measure':cached,'full_halo_Gamma_float':full,'residual_halo_Gamma_float':res,'Delta_Gamma_float':value,'margin_to_1over18_float':float(threshold)-value,'minimum_exact_cap_slack':str(min(slacks))};rec['mass_points'].append(r);counts['mass_triple_active' if active else 'mass_triple_null']+=1;top.append({'pair_index':index,'category':category,'positions_units':[Q1,Q2]}|r);top=sorted(top,key=lambda z:z['Delta_Gamma_float'],reverse=True)[:20]
   if value>float(trigger):
    triggered={'pair_index':index,'mass_point_index':len(rec['mass_points'])-1,'positions':[[str(x) for x in q1],[str(x) for x in q2]],'masses':[str(b1),str(b2)],'Delta_Gamma_float':value,'trigger':str(trigger)};break
  records.append(rec);counts[category]+=1
  if triggered:break
  if (index+1)%128==0:print(json.dumps({'event':'progress','position_pairs':index+1,'mass_points':masspoints_total,'unique_actual_measures':len(cache),'max_area_float':top[0]['Delta_Gamma_float'],'seconds':time.monotonic()-start}),flush=True)
 dp=Path(a.data);payload={'schema':'r076-two-outside-finite-screen-data-v1','parameters':params,'position_pairs':records,'actual_pair_count':len(records),'mass_point_count':masspoints_total,'triggered':triggered};dp.write_bytes(gzip.compress((json.dumps(payload,separators=(',',':'))+'\n').encode(),mtime=0))
 out={'schema':'r076-two-outside-finite-screen-summary-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'input_path':str(ip),'input_sha256':sha(ip),'parameters':params,'data_path':str(dp),'data_sha256':sha(dp),'generated_unique_position_pairs':len(pool),'actual_position_pairs_evaluated':len(records),'mass_points_checked':masspoints_total,'unique_actual_measure_area_evaluations':len(cache),'duplicate_measure_area_evaluations_avoided':masspoints_total-len(cache),'complete_exact_cap_checks':capchecks,'category_and_triple_counts':dict(counts),'triggered':triggered,'top20':top,'maximum_area_float':top[0]['Delta_Gamma_float'],'minimum_margin_float':top[0]['margin_to_1over18_float'],'seconds':time.monotonic()-start,'scope':'One purpose-designed bounded finite screen only. All fullS3 bbox caps exact, compressed to3outside-mask bounds; every sampled mass rechecked against all cap rows. All7full/3res local labels included; Gamma area floating and not certified. Triple-null diagnostic does not filter samples. Not a proof over positions/masses or a repeat search.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n')
 if triggered:
  data={'threshold':'1','background':'Lebesgue R2','geometry':'centered axis-parallel common-radius squares','S_indices':list(range(9)),'positions':[[str(F(x,C)),str(F(y,C))] for x,y in S]+triggered['positions'],'masses':['1/9']*9+triggered['masses']};Path('computation/inputs/r076/screen-trigger.json').write_text(json.dumps(data,indent=2)+'\n')
 print(json.dumps({'output':a.out,'position_pairs':len(records),'mass_points':masspoints_total,'unique_measures':len(cache),'maximum_area_float':out['maximum_area_float'],'triggered':triggered,'top':top[0],'seconds':out['seconds']},indent=2))
if __name__=='__main__':main()
