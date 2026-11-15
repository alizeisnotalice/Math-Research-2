#!/usr/bin/env python3
"""Finite exact B&B with fixed certified family point and full split ledger."""
import argparse,hashlib,heapq,json,multiprocessing,os,platform,resource,sys,time
from fractions import Fraction as F
from pathlib import Path
from mass_family_envelope import MassFamilyEnvelope
from periodic_halo import union_area

family=None
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upper(cell):
 start=time.process_time();value=F(family.upper(*cell)['upper_ratio'])
 return value,time.process_time()-start
def main():
 global family
 p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--output',required=True);p.add_argument('--reference-certificate',required=True);p.add_argument('--reference-input',required=True);p.add_argument('--reference-point',required=True,nargs=2);p.add_argument('--workers',type=int,default=6);p.add_argument('--gap',default='1/100000000');p.add_argument('--max-leaves',type=int,default=5000);p.add_argument('--seconds',type=float,default=3600);p.add_argument('--checkpoint-seconds',type=float,default=30);p.add_argument('--bits',type=int,default=40);p.add_argument('--radius',default='5');args=p.parse_args()
 goal=F(args.gap)
 if args.workers<2 or goal<=0 or args.max_leaves<1 or args.seconds<=0 or args.checkpoint_seconds<=0:p.error('invalid positive limits or fewer than two workers')
 start=time.monotonic();pcpu=time.process_time();family=MassFamilyEnvelope(json.loads(Path(args.input).read_text()),args.bits,F(args.radius));init=time.monotonic()-start;initcpu=time.process_time()-pcpu
 if not family.meta['all_corner_labels_strictly_positive']:raise ValueError('strictly positive family corner labels required')
 cert= json.loads(Path(args.reference_certificate).read_text());ref= json.loads(Path(args.reference_input).read_text());point=family.point(*args.reference_point)
 if cert['input_sha256']!=sha(args.reference_input) or cert['code_sha256']!=sha(Path(__file__).with_name('periodic_halo.py')):raise ValueError('reference certificate mismatch')
 if any(F(ref[k])!=F(point[k]) for k in ('Px','Py','lam')) or len(ref['pts'])!=len(point['pts']) or any(any(F(x)!=F(y) for x,y in zip(u,v)) for u,v in zip(ref['pts'],point['pts'])):raise ValueError('reference point is not claimed family input')
 lower=F(cert['inner_ratio']);assert F(cert['mass'])==family.mass
 root=(family.a0,family.a1,family.b0,family.b1);root_upper,root_cpu=upper(root);assert root_upper>=lower
 nodes=[{'id':0,'parent':None,'cell':root,'upper':root_upper,'status':'active','children':[]}];heap=[(-root_upper,0)];pruned=[];worker_cpu=F(0);splits=0;reason='time_limit';next_save=time.monotonic()+args.checkpoint_seconds;before=resource.getrusage(resource.RUSAGE_CHILDREN)
 def save(status):
  leaves=[n for n in nodes if n['status']!='split'];rects=[n['cell'] for n in leaves];domain=(family.a1-family.a0)*(family.b1-family.b0)
  assert all(family.a0<=a<b<=family.a1 and family.b0<=c<d<=family.b1 for a,b,c,d in rects)
  assert sum(((b-a)*(d-c) for a,b,c,d in rects),F(0))==union_area(rects)==domain
  assert all(n['upper']<lower for n in leaves if n['status']=='pruned')
  maxupper=max([lower]+[n['upper'] for n in leaves]);difference=maxupper-lower
  if status=='target_gap':assert difference<=goal
  encoded=[{'id':n['id'],'parent':n['parent'],'cell':list(map(str,n['cell'])),'upper':str(n['upper']),'status':n['status'],'children':n['children']} for n in nodes]
  result={'status':status,'scope':'pending independent review; fixed certified point lower, per-box corner mass upper','input_sha256':sha(args.input),'driver_sha256':sha(__file__),'envelope_sha256':sha(Path(__file__).with_name('mass_family_envelope.py')),'reference_input_sha256':sha(args.reference_input),'reference_certificate_sha256':sha(args.reference_certificate),'reference_point':args.reference_point,'python':sys.version,'platform':platform.platform(),'family':family.meta,'workers':args.workers,'target_gap':str(goal),'max_leaves':args.max_leaves,'initialization_seconds':init,'initialization_cpu_seconds':initcpu,'elapsed_seconds':time.monotonic()-start,'parent_cpu_seconds':time.process_time()-pcpu,'aggregate_worker_cpu_seconds':float(worker_cpu)+root_cpu,'splits':splits,'leaf_count':len(leaves),'partition_area_and_union_verified':True,'lower':str(lower),'upper':str(maxupper),'gap':str(difference),'lower_float':float(lower),'upper_float':float(maxupper),'gap_float':float(difference),'nodes':encoded}
  temp=Path(args.output+'.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');os.replace(temp,args.output)
  print(json.dumps({k:v for k,v in result.items() if k!='nodes'}),flush=True)
 pool=multiprocessing.get_context('fork').Pool(args.workers)
 try:
  while heap:
   if -heap[0][0]-lower<=goal:reason='target_gap';break
   nleaves=len(heap)+len(pruned)
   if nleaves>=args.max_leaves:reason='leaf_limit';break
   if time.monotonic()-start>=args.seconds:break
   count=min(max(1,args.workers//2),len(heap),args.max_leaves-nleaves);selected=[heapq.heappop(heap)[1] for _ in range(count)];children=[]
   for ident in selected:
    a,b,c,d=nodes[ident]['cell'];am,bm=(a+b)/2,(c+d)/2
    if (b-a)/(family.a1-family.a0)>=(d-c)/(family.b1-family.b0):children.extend(((a,am,c,d),(am,b,c,d)))
    else:children.extend(((a,b,c,bm),(a,b,bm,d)))
   values=pool.map(upper,children,chunksize=1)
   for pos,ident in enumerate(selected):
    parent=nodes[ident];parent['status']='split';kids=[]
    for cell,(value,usedcpu) in zip(children[2*pos:2*pos+2],values[2*pos:2*pos+2]):
     worker_cpu+=F.from_float(usedcpu);assert value<=parent['upper'];idx=len(nodes);state='pruned' if value<lower else 'active';nodes.append({'id':idx,'parent':ident,'cell':cell,'upper':value,'status':state,'children':[]});kids.append(idx)
     if state=='pruned':pruned.append(idx)
     else:heapq.heappush(heap,(-value,idx))
    parent['children']=kids;splits+=1
   if time.monotonic()>=next_save:save('running');next_save=time.monotonic()+args.checkpoint_seconds
  else:reason='target_gap'
 finally:pool.close();pool.join()
 save(reason)
if __name__=='__main__':main()
