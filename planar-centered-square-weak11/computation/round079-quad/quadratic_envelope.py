#!/usr/bin/env python3
"""Exact cell certificate for two specified continuous quadratic envelopes.

This estimates E_quad minus the actual residual halo. E_quad is an OUTER
state for old-involving actual halo witnesses; its inner area is NOT an
actual halo lower bound. No parameter search is implemented.
"""
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, gzip, hashlib, heapq, json, math, os, platform, sys, time

HALF=F(1,2)
EXPECTED={
 'corrected':'d4a87c2ef5f82165467bed876c71a0a86d19dce89012ea37b76bf83583fdd15a',
 'farface':'e343e448b9582a731ddaccca37452adcce1c6f28b4153aa520a0ec2b4a3ef574',
}
VARIANTS=('exact_pureold_plus_three_mixed','four_soft_branches')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rect_area(r): return (r[1]-r[0])*(r[3]-r[2]) if r else F(0)
def intersect(*rects):
 if any(r is None for r in rects): return None
 a=max(r[0] for r in rects);b=min(r[1] for r in rects)
 c=max(r[2] for r in rects);d=min(r[3] for r in rects)
 return (a,b,c,d) if a<b and c<d else None
def area_minus(rect,residual):
 """Exact inclusion-exclusion for a rectangle minus at most three rectangles."""
 if rect is None:return F(0)
 rr=[z for r in residual if (z:=intersect(rect,r)) is not None]
 ans=rect_area(rect)
 for count in range(1,len(rr)+1):
  for part in combinations(rr,count):ans+=(-1 if count%2 else 1)*rect_area(intersect(*part))
 assert 0<=ans<=rect_area(rect)
 return ans
def source_bbox(points):
 return (min(p[0] for p in points),max(p[0] for p in points),min(p[1] for p in points),max(p[1] for p in points))
def radius_range(cell,bbox):
 """Min/max over a positive-coordinate cell of farthest source distance."""
 u,v,s,t=cell;a,b,c,d=bbox
 lo=max((b-a)/2,(d-c)/2,b-v,u-a,d-t,s-c,F(0))
 hi=max(b-u,v-a,d-s,t-c,F(0))
 assert 0<=lo<=hi
 return lo,hi
def qmin(X,Y,rlo,n):
 """Closed-interval minimum on [rlo,1/2]; all quantities rational."""
 assert rlo<=HALF
 c=F(n*n,(n+1)*(n+1));h=F(n+1,n*n)
 A=4-c;B=c*(X+Y-2*h);C=-c*(X-h)*(Y-h)
 r0=-B/(2*A);star=max(rlo,min(HALF,r0));value=A*star*star+B*star+C
 return value,star
def branch_state(cell,bbox,B,n):
 """in/out/unknown refer only to the continuous quadratic outer state."""
 rlo,rhi=radius_range(cell,bbox)
 if rlo>=HALF:return 'out'
 low,_=qmin(cell[0],cell[2],rlo,n)
 if low>=B:return 'out'
 if rhi<HALF:
  high,_=qmin(cell[1],cell[3],rhi,n)
  if high<B:return 'in'
 return 'unknown'
def union_state(states):
 if 'in' in states:return 'in'
 return 'out' if all(x=='out' for x in states) else 'unknown'
def split(cell):
 u,v,s,t=cell;x=(u+v)/2;y=(s+t)/2
 return [(u,x,s,y),(x,v,s,y),(u,x,y,t),(x,v,y,t)]
def make_residual(points,masses,bits,n):
 """Three genuine outside labels; strict receivers, closed area enclosures."""
 den=math.lcm(n*n,*[x.denominator for p in points for x in p],*[m.denominator for m in masses]);scale=den*(1<<bits)
 inner=[];outer=[];catalog=[]
 for mask in [1,2,3]:
  jj=[j for j in range(2) if mask>>j&1];mass=sum(masses[j] for j in jj);box=source_bbox([points[j] for j in jj]);span=max(box[1]-box[0],box[3]-box[2]);row={'outside_mask':mask,'mass':str(mass),'bbox':list(map(str,box)),'eligible':mass>span*span}
  if row['eligible']:
   z=mass*scale*scale;assert z.denominator==1
   k=math.isqrt(z.numerator//4);up=k if 4*k*k==z else k+1;assert 4*k*k<=z<=4*up*up
   a,b,c,d=box;ri=(b-F(k,scale),a+F(k,scale),d-F(k,scale),c+F(k,scale));ro=(b-F(up,scale),a+F(up,scale),d-F(up,scale),c+F(up,scale))
   assert ri[0]<ri[1] and ri[2]<ri[3]
   inner.append(ri);outer.append(ro);row.update(radius_lower=str(F(k,scale)),radius_upper=str(F(up,scale)),inner=list(map(str,ri)),outer=list(map(str,ro)))
  catalog.append(row)
 return inner,outer,catalog
def validate_input(path,case_id,bits):
 path=Path(path);assert sha(path)==EXPECTED[case_id]
 data=json.loads(path.read_text());p=[tuple(map(F,x)) for x in data['positions']];mass=list(map(F,data['masses']));n=20
 assert F(data['threshold'])==1 and p[:400]==[(F(21*i,400),F(21*j,400)) for i in range(20) for j in range(20)] and mass[:400]==[F(1,400)]*400
 assert data['S_indices']==list(range(400)) and len(p)==402
 if case_id=='corrected':assert p[400:]==[(F(-71,200),F(0)),(F(0),F(-171,200))] and mass[400:]==[F(1,4),F(1,2)]
 else:assert p[400:]==[(F(-1,2),F(0)),(F(0),F(-1,2))] and mass[400:]==[F(37,100)]*2
 # Reflect source and receiver coordinates together: X=-x,Y=-y.
 positive=[(-x,-y) for x,y in p[400:]];origin=(F(0),F(0));branches=[]
 for mask in [0,1,2,3]:
  jj=[j for j in range(2) if mask>>j&1]
  branches.append((source_bbox([origin]+[positive[j] for j in jj]),sum((mass[400+j] for j in jj),F(0))))
 ri,ro,cat=make_residual(positive,mass[400:],bits,n)
 return {'case_id':case_id,'input_path':str(path),'input_sha256':sha(path),'n':n,'branches':branches,'residual_inner':ri,'residual_outer':ro,'residual_catalog':cat,'pureold':(F(1,400),F(1,40),F(1,400),F(1,40))}
def classify(cell,ctx,path):
 ai=area_minus(cell,ctx['residual_inner']);ao=area_minus(cell,ctx['residual_outer']);assert ao<=ai
 # No D-state area remains when an actual residual inner rectangle union
 # already covers this cell. Avoid all four polynomial tests in that case.
 if ai==0:
  return {'path':path,'branch_states':['skipped_residual_cover']*4,'mixed_state':'residual_covered','soft_state':'residual_covered','remaining_residual_inner':ai,'remaining_residual_outer':ao,'known_pureold_remaining_lower':F(0),'known_pureold_remaining_upper':F(0),'lower':[F(0),F(0)],'upper':[F(0),F(0)],'unknown_area':[F(0),F(0)]}
 states=[branch_state(cell,bbox,B,ctx['n']) for bbox,B in ctx['branches']];mixed=union_state(states[1:]);soft=union_state(states)
 old=intersect(cell,ctx['pureold']);oldlo=area_minus(old,ctx['residual_outer']);oldhi=area_minus(old,ctx['residual_inner'])
 if mixed=='in':exactlo,exacthi=ao,ai
 elif mixed=='out':exactlo,exacthi=oldlo,oldhi
 else:exactlo,exacthi=oldlo,ai
 softlo=ao if soft=='in' else F(0);softhi=F(0) if soft=='out' else ai
 lower=[exactlo,softlo];upper=[exacthi,softhi];assert all(0<=l<=u<=rect_area(cell) for l,u in zip(lower,upper))
 return {'path':path,'branch_states':states,'mixed_state':mixed,'soft_state':soft,'remaining_residual_inner':ai,'remaining_residual_outer':ao,'known_pureold_remaining_lower':oldlo,'known_pureold_remaining_upper':oldhi,'lower':lower,'upper':upper,'unknown_area':[rect_area(cell) if mixed=='unknown' else F(0),rect_area(cell) if soft=='unknown' else F(0)]}
def serial_leaf(leaf):
 def conv(v):
  if isinstance(v,F):return str(v)
  if isinstance(v,list):return [conv(x) for x in v]
  return v
 return {k:conv(v) for k,v in leaf.items()}
def tile_run(task):
 """One disjoint initial tile, no nested workers or multithreaded BLAS."""
 ctx,cell,tile_id,tol,maxnodes,maxdepth=task;t=time.monotonic();root=classify(cell,ctx,'');leaves={'':root};boxes={'':cell};heap=[];created=1;split_count=0;depth=0
 lower=root['lower'][:];upper=root['upper'][:]
 def enqueue(leaf):
  gap=max(u-l for l,u in zip(leaf['lower'],leaf['upper']))
  if gap>0:heapq.heappush(heap,(-gap,leaf['path']))
 enqueue(root);reason='target_met'
 while max(u-l for l,u in zip(lower,upper))>tol:
  if not heap:reason='no_refinable_leaf';break
  _,path=heapq.heappop(heap);leaf=leaves[path]
  if len(path)>=maxdepth:reason='maximum_depth';break
  if created+4>maxnodes:reason='maximum_nodes';break
  box=boxes.pop(path);del leaves[path]
  for i in range(2):lower[i]-=leaf['lower'][i];upper[i]-=leaf['upper'][i]
  for j,sub in enumerate(split(box)):
   key=path+str(j);child=classify(sub,ctx,key);leaves[key]=child;boxes[key]=sub;enqueue(child)
   for i in range(2):lower[i]+=child['lower'][i];upper[i]+=child['upper'][i]
  created+=4;split_count+=1;depth=max(depth,len(path)+1)
 # Complete leaves, including every unresolved branch, remain in output.
 assert len(leaves)==1+3*split_count and sum((rect_area(b) for b in boxes.values()),F(0))==rect_area(cell)
 assert lower==[sum((z['lower'][i] for z in leaves.values()),F(0)) for i in range(2)] and upper==[sum((z['upper'][i] for z in leaves.values()),F(0)) for i in range(2)]
 return {'case_id':ctx['case_id'],'tile_id':tile_id,'root_box':list(map(str,cell)),'pid':os.getpid(),'termination':reason,'target_per_variant':str(tol),'created_node_count':created,'split_count':split_count,'leaf_count':len(leaves),'maximum_depth':depth,'lower':list(map(str,lower)),'upper':list(map(str,upper)),'gaps':[str(u-l) for l,u in zip(lower,upper)],'unresolved_leaf_area':[str(sum((z['unknown_area'][i] for z in leaves.values()),F(0))) for i in range(2)],'leaf_records':[serial_leaf(leaves[key]) for key in sorted(leaves)],'seconds':time.monotonic()-t}
def encode_gzip(path,data):
 with Path(path).open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as z:z.write((json.dumps(data,separators=(',',':'))+'\n').encode())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--parameters',required=True);ap.add_argument('--output-dir',required=True);ap.add_argument('--workers',type=int,default=1);args=ap.parse_args();start=time.monotonic()
 assert Path('/Volumes/data').is_mount(),'data not mounted'
 pp=Path(args.parameters).resolve();params=json.loads(pp.read_text());assert params['schema']=='r079-two-fixed-quadratic-envelopes-v1';out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
 assert not (out/'summary.json').exists(),'refusing to replace previous production outputs'
 assert 1<=args.workers<=12
 cases=params['cases'];assert {c['id'] for c in cases}==set(EXPECTED) and len(cases)==2
 bits=int(params['residual_bits']);assert bits>=72;tiles=int(params['tiles_per_axis']);assert 1<=tiles<=4
 tolerance=F(params['target_area_gap']);assert tolerance>0
 e=F(1,400);axis=[e+(HALF-e)*i/tiles for i in range(tiles+1)];tasks=[];contexts={}
 for case in cases:
  ctx=validate_input(pp.parent/case['input'],case['id'],bits);contexts[case['id']]=ctx
  for i in range(tiles):
   for j in range(tiles):tasks.append((ctx,(axis[i],axis[i+1],axis[j],axis[j+1]),f'{i}-{j}',tolerance/(tiles*tiles),int(params['maximum_nodes_per_tile']),int(params['maximum_depth'])))
 if args.workers==1:results=[tile_run(task) for task in tasks]
 else:
  with ProcessPoolExecutor(max_workers=args.workers) as pool:results=list(pool.map(tile_run,tasks,chunksize=1))
 summaries=[]
 for case in cases:
  cid=case['id'];rr=sorted((z for z in results if z['case_id']==cid),key=lambda z:z['tile_id']);ctx=contexts[cid]
  lower=[sum((F(z['lower'][i]) for z in rr),F(0)) for i in range(2)];upper=[sum((F(z['upper'][i]) for z in rr),F(0)) for i in range(2)];gaps=[u-l for l,u in zip(lower,upper)]
  cert={'schema':'r079-quadratic-cell-certificate-v1','case_id':cid,'input_sha256':ctx['input_sha256'],'source_sha256':sha(__file__),'parameters_sha256':sha(pp),'positive_Gamma':['1/400','1/2','1/400','1/2'],'variants':list(VARIANTS),'n':20,'h':'21/400','c':'400/441','residual_bits':bits,'residual_catalog':ctx['residual_catalog'],'pureold_rectangle':list(map(str,ctx['pureold'])),'children_order':['lowerX_lowerY','upperX_lowerY','lowerX_upperY','upperX_upperY'],'tiles':rr,'scope':'All in/inner statements refer to continuous quadratic envelope state only. This is not an actual halo lower certificate.'}
  cp=out/f'{cid}-certificate.json.gz';assert not cp.exists();encode_gzip(cp,cert)
  summaries.append({'case_id':cid,'input_sha256':ctx['input_sha256'],'certificate':cp.name,'certificate_sha256':sha(cp),'certificate_bytes':cp.stat().st_size,'variants':{name:{'envelope_minus_true_residual_lower':str(lower[i]),'envelope_minus_true_residual_upper':str(upper[i]),'lower_decimal':float(lower[i]),'upper_decimal':float(upper[i]),'gap':str(gaps[i]),'target_met':gaps[i]<=tolerance,'unresolved_leaf_area':str(sum((F(z['unresolved_leaf_area'][i]) for z in rr),F(0))),'actual_deletion_upper_valid':True,'this_lower_is_not_actual_deletion_lower':True} for i,name in enumerate(VARIANTS)},'tile_count':len(rr),'leaf_count':sum(z['leaf_count'] for z in rr),'created_node_count':sum(z['created_node_count'] for z in rr),'all_tiles_target_met':all(z['termination']=='target_met' for z in rr),'worker_pids':sorted({z['pid'] for z in rr})})
 summary={'schema':'r079-quadratic-envelope-summary-v1','source_sha256':sha(__file__),'parameters_sha256':sha(pp),'pid':os.getpid(),'python':sys.version,'platform':platform.platform(),'workers':args.workers,'worker_pids':sorted({z['pid'] for z in results}),'target_area_gap':str(tolerance),'cases':summaries,'seconds':time.monotonic()-start,'scope':'Two fixed actual inputs, continuous quadratic oldcapture OUTER states minus identical complete actual residual. State-inner lower bounds cannot be used as actual halo lower bounds. No parameter search or general C2 claim.'}
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
