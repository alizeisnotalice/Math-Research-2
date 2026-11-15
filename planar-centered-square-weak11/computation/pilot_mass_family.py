#!/usr/bin/env python3
"""Uniform-grid upper-bound pilot; parallel workers share one prepared family."""
import argparse,hashlib,json,multiprocessing,platform,resource,statistics,sys,time
from fractions import Fraction as F
from pathlib import Path
from mass_family_envelope import MassFamilyEnvelope

family=None
def one(cell):
 start=time.monotonic();cpu=time.process_time();result=family.upper(*cell)
 return {'cell':list(map(str,cell)),'upper':result['upper_ratio'],'seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu}
def main():
 global family
 p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--output',required=True);p.add_argument('--reference-lower');p.add_argument('--reference-certificate');p.add_argument('--reference-input');p.add_argument('--reference-point',nargs=2);p.add_argument('--workers',type=int,default=1);p.add_argument('--bits',type=int,default=40);p.add_argument('--radius',default='5');p.add_argument('--allow-zero-corners',action='store_true',help='valid upper bounds only; no strict-positivity convergence budget');args=p.parse_args()
 if args.workers<1:p.error('workers must be positive')
 if bool(args.reference_lower)==bool(args.reference_certificate):p.error('choose exactly one reference source')
 if args.reference_certificate and (not args.reference_input or not args.reference_point):p.error('certificate requires reference input and point')
 start=time.monotonic();cpu_start=time.process_time();family=MassFamilyEnvelope(json.loads(Path(args.input).read_text()),args.bits,F(args.radius));init=time.monotonic()-start;init_cpu=time.process_time()-cpu_start;grids=[];last=None;binding=None
 if not family.meta['all_corner_labels_strictly_positive'] and not args.allow_zero_corners:raise ValueError('strictly positive corner labels required unless explicitly overridden')
 if args.reference_certificate:
  certpath,refpath=Path(args.reference_certificate),Path(args.reference_input);cert=json.loads(certpath.read_text());refdata=json.loads(refpath.read_text());candidate=family.point(*args.reference_point)
  refsha=hashlib.sha256(refpath.read_bytes()).hexdigest()
  if cert['input_sha256']!=refsha or cert['code_sha256']!=hashlib.sha256(Path(__file__).with_name('periodic_halo.py').read_bytes()).hexdigest():raise ValueError('reference certificate input or code mismatch')
  if any(F(refdata[k])!=F(candidate[k]) for k in ('Px','Py','lam')) or len(refdata['pts'])!=len(candidate['pts']) or any(any(F(x)!=F(y) for x,y in zip(u,v)) for u,v in zip(refdata['pts'],candidate['pts'])):raise ValueError('reference input differs from claimed family point')
  ref=F(cert['inner_ratio']);binding={'input_sha256':refsha,'certificate_sha256':hashlib.sha256(certpath.read_bytes()).hexdigest(),'point':args.reference_point}
 else:ref=F(args.reference_lower)
 before=resource.getrusage(resource.RUSAGE_CHILDREN)
 pool=multiprocessing.get_context('fork').Pool(args.workers) if args.workers>1 else None
 try:
  for n in (1,2,4,8):
   cells=[(family.a0+(family.a1-family.a0)*i/n,family.a0+(family.a1-family.a0)*(i+1)/n,family.b0+(family.b1-family.b0)*j/n,family.b0+(family.b1-family.b0)*(j+1)/n) for i in range(n) for j in range(n)]
   tick=time.monotonic();results=pool.map(one,cells,chunksize=1) if pool else list(map(one,cells));elapsed=time.monotonic()-tick;upper=max(F(x['upper']) for x in results)
   assert last is None or upper<=last
   gaps=upper-ref;durations=[v['seconds'] for v in results]
   grids.append({'partitions':n*n,'upper':str(upper),'upper_float':float(upper),'gap_to_reference':str(gaps),'gap_float':float(gaps),'wall_seconds':elapsed,'mean_cell_seconds':statistics.mean(durations),'median_cell_seconds':statistics.median(durations),'min_cell_seconds':min(durations),'max_cell_seconds':max(durations),'cells':results});last=upper
 finally:
  if pool:pool.close();pool.join()
 after=resource.getrusage(resource.RUSAGE_CHILDREN)
 result={'input_sha256':hashlib.sha256(Path(args.input).read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'envelope_sha256':hashlib.sha256(Path(__file__).with_name('mass_family_envelope.py').read_bytes()).hexdigest(),'python':sys.version,'platform':platform.platform(),'workers':args.workers,'family':family.meta,'reference_lower':str(ref),'reference_binding':binding,'initialization_seconds':init,'initialization_cpu_seconds':init_cpu,'parent_cpu_seconds':time.process_time()-cpu_start,'children_user_cpu_seconds':after.ru_utime-before.ru_utime,'children_system_cpu_seconds':after.ru_stime-before.ru_stime,'elapsed_seconds':time.monotonic()-start,'grids':grids}
 Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='grids'}|{'grids':[{k:v for k,v in g.items() if k!='cells'} for g in grids]},indent=2))
if __name__=='__main__':main()
