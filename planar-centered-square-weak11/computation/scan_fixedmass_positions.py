#!/usr/bin/env python3
"""R8-based labelled geometry A/B grid with an independently fixed mass vector."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from scan_supercell_2x2 import canon,digest
from scan_supercell_joint import make
from scan_supercell_mass import target_delta
from periodic_halo_dominated import calculate_dominated

def worker(t,base):
 start=time.monotonic();R=5
 try:
  while True:
   try:r=calculate_dominated(t['input'],base,[2,2],F(t['delta']),40,F(R),require_mass_conservation=True);break
   except ValueError as e:
    if str(e)!='transport tail bound does not prove the requested radius cutoff':raise
    R+=1
  return {'id':t['id'],'input_sha256':t['input_sha256'],'delta':t['delta'],'kappa':r['transport_matching']['kappa'],'tail_radius':R,'tail_bound':r['tail_bound'],'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates']}
 except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--base',required=True,type=Path);ap.add_argument('--mass-input',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--workers',type=int,default=6);args=ap.parse_args()
 if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
 base=json.loads(args.base.read_text());massdata=json.loads(args.mass_input.read_text());n=len(base['pts'])
 if len(massdata['pts'])!=4*n:raise ValueError('wrong labelled input length')
 if sum(F(p[2]) for p in massdata['pts'])!=4*sum(F(p[2]) for p in base['pts']):raise ValueError('mass not conserved')
 masses=[F(p[2]) for p in massdata['pts']];unique={}
 for i,j in product(range(-8,9),repeat=2):
  a=F(19,16)+F(i,32);b=F(27,64)+F(j,32);d,_=make(base,(a,b,F(0)));d['pts']=[[p[0],p[1],str(m)] for p,m in zip(d['pts'],masses)];delta=target_delta(base,d);h=digest(d);alias={'i':i,'j':j,'a':str(a),'b':str(b)}
  if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'delta':str(delta),'aliases':[]}
  elif canon(unique[h]['input'])!=canon(d) or unique[h]['delta']!=str(delta):raise RuntimeError('collision')
  unique[h]['aliases'].append(alias)
 tasks=list(unique.values());args.out.mkdir(parents=True,exist_ok=True);manifest={'schema':'fixed-mass-two-position-modes-v1','base_sha256':sha(args.base),'mass_input_sha256':sha(args.mass_input),'dominated_kernel_sha256':sha(Path(__file__).with_name('periodic_halo_dominated.py')),'geometry_builder_sha256':sha(Path(__file__).with_name('scan_supercell_joint.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'grid':'A=19/16+i/32, B=27/64+j/32, i,j=-8..8','bits':40,'initial_tail_radius':5,'labels':'canonical (copy_x,copy_y,base_index)','workers':args.workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks};(args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'pid={os.getpid()} points={len(tasks)} workers={args.workers}',flush=True)
 start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
 with (args.out/'completed.jsonl').open('x') as stream,cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
  futs={pool.submit(worker,t,base):t for t in tasks}
  for f in cf.as_completed(futs):
   try:r=f.result()
   except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
   results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
   if now-last>=60:print(f'completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
 results.sort(key=lambda r:r['id']);(args.out/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
 for r in top:
  t=byid[r['id']];r['aliases']=t['aliases'];(args.out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
 u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'points':len(tasks),'completed':len(results)-len(errors),'errors':len(errors),'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(args.out/'results.json'),'top':top};(args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'done={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
 if errors:raise RuntimeError('candidate failure')
if __name__=='__main__':main()
