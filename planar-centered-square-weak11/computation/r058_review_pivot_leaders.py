#!/usr/bin/env python3
"""Standalone rational-eligibility replay of three maximizers, no C++ imports."""
from pathlib import Path
import hashlib,itertools,json,math,os,platform,sys,time
ROOT=Path(__file__).resolve().parents[1];DIR=ROOT/'computation/results/r058-pivot-search';B=3.299999995647862

def union(rect):
 # Vertical rather than horizontal sweep; half-open boundaries immaterial.
 ys=sorted(set(t for r in rect for t in r[2:]));ans=0
 for a,b in zip(ys,ys[1:]):
  ints=sorted((x0,x1)for x0,x1,y0,y1 in rect if y0<=a and b<=y1 and x0<x1);stop=-math.inf;length=0
  for l,h in ints:
   if h>stop:length+=h-max(l,stop);stop=h
  ans+=(b-a)*length
 return ans

def replay(rec):
 c=rec['input'];pts=c['points_xy_integer_weights'];n=c['n'];Q=c['coordinate_denominator'];W=sum(p[2]for p in pts);full=(1<<n)-1;stats={}
 for s in range(1,1<<n):
  ids=[j for j in range(n)if s>>j&1];x=[pts[j][0]for j in ids];y=[pts[j][1]for j in ids];m=sum(pts[j][2]for j in ids);d=max(max(x)-min(x),max(y)-min(y));stats[s]=(m,min(x),max(x),min(y),max(y),d,len(ids))
 eligible={s for s,t in stats.items()if t[5]*t[5]*W<t[0]*Q*Q}
 cacheT={}
 def cap(j,s):
  key=j,s
  if key not in cacheT:cacheT[key]=max(stats[t][0]for t in eligible if t&1<<j and t&s==t)
  return cacheT[key]
 def square(j,m):
  h=math.sqrt(m/W)/2;x,y,_=pts[j];return (x/Q-h,x/Q+h,y/Q-h,y/Q+h)
 def capture(s):
  m,x0,x1,y0,y1,_,_=stats[s];h=math.sqrt(m/W)/2;return (x1/Q-h,x0/Q+h,y1/Q-h,y0/Q+h)
 leafcache={}
 def leaf(s):
  if s in leafcache:return leafcache[s]
  ids=[j for j in range(n)if s>>j&1];m,x0,x1,y0,y1,d,_=stats[s];rank=max(stats[t][6]for t in eligible if t&s==t);u=min(B*max(m/W,d*d/Q/Q),rank*m/W);g=[square(j,pts[j][2])for j in ids]
  if s in eligible:g.append(capture(s))
  for key in [lambda j:(pts[j][2],j),lambda j:(-pts[j][2],j),lambda j:(pts[j][0],j),lambda j:(pts[j][1],j)]:
   prefix=0;r=[]
   for j in sorted(ids,key=key):prefix|=1<<j;r.append(square(j,cap(j,prefix)))
   u=min(u,union(r))
  leafcache[s]=(u,g,union(g));return leafcache[s]
 best=math.inf;bd=None;bl=None;checked=0
 for k in range(4):
  for deleted in itertools.combinations(range(n),k):
   D=sum(1<<j for j in deleted);rem=full^D;edges=[s for s in eligible if s&rem==s];parts=[]
   while rem:
    component=rem&-rem
    while True:
     enlarged=component
     for s in edges:
      if s&component:enlarged|=s
     if enlarged==component:break
     component=enlarged
    parts.append(component);rem^=component
   rr=[square(j,cap(j,full))for j in deleted];resid=0
   for s in parts:
    u,g,ga=leaf(s);rr.extend(g);resid+=u-ga
   z=union(rr)+resid;checked+=1
   if z<best:best=z;bd=D;bl=parts
 halo=union([capture(s)for s in eligible]);assert abs(best-rec['score'])<1e-10;assert best+1e-10>=halo
 return {'input':c,'reference_score':rec['score'],'replayed_best_score':best,'best_deleted_mask':bd,'best_leaf_masks':bl,'deletions_checked':checked,'actual_all_subset_halo_float':halo,'overestimate':best-halo,'eligible_count':len(eligible),'full_rank':max(stats[s][6]for s in eligible),'bbox_d_squared':stats[full][5]**2/Q**2,'input_sha256':hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()}

def main():
 assert os.path.ismount('/Volumes/data');start=time.perf_counter();sources=[];results=[]
 for n in (8,10,12):
  pairs=[]
  for p in DIR.glob(f'n{n}-s*-k3.json'):
   d=json.loads(p.read_text())
   if d['count']==2000:pairs.append((d['leaders'][0]['score'],p,d['leaders'][0]))
  _,p,r=max(pairs,key=lambda q:q[0]);sources.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});results.append(replay(r));print('checked',n,results[-1]['replayed_best_score'],results[-1]['actual_all_subset_halo_float'],flush=True)
 result={'status':'PASS three complete deletion and actual halo replays','scope':'Floating rectangle areas with independently exact integer eligibility; numerical evidence only, not universal bound.','results':results,'sources':sources,'python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'elapsed_seconds':time.perf_counter()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};p=DIR/'independent-three-leader-replay.json';p.write_text(json.dumps(result,indent=2)+'\n');print('result_sha256',hashlib.sha256(p.read_bytes()).hexdigest());print('elapsed',result['elapsed_seconds'])
if __name__=='__main__':main()
