#!/usr/bin/env python3
"""Single-quadrant original-coordinate-group moves from fixed R25 target."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
from periodic_halo_dominated import calculate_dominated
from scan_supercell_2x2 import canon,digest
from scan_supercell_mass import target_delta

def moved(base,target,quad,axis,members,shift):
 d={k:v for k,v in target.items()};pts=[[F(v) for v in p] for p in target['pts']];n=len(base['pts']);P=F(base['Px' if axis==0 else 'Py']);superP=F(target['Px' if axis==0 else 'Py'])
 for i in members:pts[quad*n+i][axis]=(pts[quad*n+i][axis]+F(shift)*P)%superP
 d['pts']=[[str(v) for v in p] for p in pts];return d

def worker(t,base):
 start=time.monotonic();R=5
 try:
  while True:
   try:r=calculate_dominated(t['input'],base,[2,2],F(t['delta']),40,F(R),require_mass_conservation=True);break
   except ValueError as e:
    if str(e)!='transport tail bound does not prove the requested radius cutoff':raise
    R+=1
  return {'id':t['id'],'input_sha256':t['input_sha256'],'delta':t['delta'],'kappa':r['transport_matching']['kappa'],'tail_radius':R,'tail_bound':r['tail_bound'],'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates'],'window_atoms_distinct':r['window_atoms_distinct'],'box_counts':r['counts'],'inner_rectangles_unique':r['inner_rectangles_unique'],'outer_rectangles_unique':r['outer_rectangles_unique']}
 except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);ap.add_argument('--target',type=Path,required=True);ap.add_argument('--previous',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=6);args=ap.parse_args()
 if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
 base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());n=len(base['pts']);kernel=sha(Path(__file__).with_name('periodic_halo_dominated.py'));pm=json.loads((args.previous/'manifest.json').read_text());pr=json.loads((args.previous/'results.json').read_text())
 if pm['base_sha256']!=sha(args.base) or pm['dominated_kernel_sha256']!=kernel:raise ValueError('prior incompatible')
 unique={};cache={};rr={r['id']:r for r in pr if 'error' not in r}
 for t in pm['tasks']:
  r=rr.get(t['id'])
  if r and t['input_sha256']==r['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r)
 def add(d,alias):
  h=digest(d);delta=target_delta(base,d)
  if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'delta':str(delta),'aliases':[]}
  elif canon(unique[h]['input'])!=canon(d) or unique[h]['delta']!=str(delta):raise RuntimeError('collision')
  unique[h]['aliases'].append(alias)
 add(target,{'kind':'baseline'});groups=[]
 for axis in (0,1):
  gs={}
  for i,p in enumerate(base['pts']):gs.setdefault(F(p[axis]),[]).append(i)
  groups.extend({'axis':axis,'coordinate':str(v),'members':members} for v,members in sorted(gs.items()))
 eps=['1/20000','1/100000','1/500000'] if 1+4*len(groups)*6<=900 else ['1/100000','1/500000']
 for quadrant in range(4):
  for g in groups:
   for e in eps:
    for sign in (-1,1):
     shift=sign*F(e);d=moved(base,target,quadrant,g['axis'],g['members'],shift);add(d,{'kind':'quadrant_group_shift','quadrant':quadrant,'group':g,'shift':str(shift)})
 tasks=list(unique.values());reuse={}
 for t in tasks:
  if t['input_sha256'] in cache:
   old,r=cache[t['input_sha256']]
   if canon(old['input'])!=canon(t['input']) or old['delta']!=t['delta']:raise ValueError('cache mismatch')
   reuse[t['id']]={**r,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':sha(args.previous/'results.json')}
 args.out.mkdir(parents=True,exist_ok=True);manifest={'schema':'quadrant-group-geometry-v1','base_sha256':sha(args.base),'target_sha256':sha(args.target),'source_prior_manifest_sha256':sha(args.previous/'manifest.json'),'source_prior_results_sha256':sha(args.previous/'results.json'),'dominated_kernel_sha256':kernel,'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'groups':groups,'eps':eps,'bits':40,'initial_radius':5,'labels':'canonical (copy_x,copy_y,base_index)','workers':args.workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks};(args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'pid={os.getpid()} points={len(tasks)} reused={len(reuse)} groups={len(groups)} eps={eps}',flush=True)
 start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
 with (args.out/'completed.jsonl').open('x') as stream,cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
  for r in reuse.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
  futs={pool.submit(worker,t,base):t for t in tasks if t['id'] not in reuse}
  for f in cf.as_completed(futs):
   try:r=f.result();r['reused']=False
   except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
   results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
   if now-last>=60:print(f'completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
 results.sort(key=lambda r:r['id']);(args.out/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
 for r in top:
  t=byid[r['id']];r['aliases']=t['aliases'];(args.out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
 # Earlier geometry scanner lacks box counts; calculate same input once to provide comparison.
 baseline=worker(tasks[0],base)
 if 'error' in baseline:raise RuntimeError('baseline diagnostic failure')
 diagnostics=[]
 for r in results:
  if 'box_counts' not in r:continue
  diagnostics.append({'id':r['id'],'box_count_changes':{k:r['box_counts'][k]-baseline['box_counts'][k] for k in r['box_counts']},'inner_rectangles_change':r['inner_rectangles_unique']-baseline['inner_rectangles_unique'],'outer_rectangles_change':r['outer_rectangles_unique']-baseline['outer_rectangles_unique'],'window_x_change':r['window_x_coordinates']-baseline['window_x_coordinates'],'window_y_change':r['window_y_coordinates']-baseline['window_y_coordinates']})
 (args.out/'geometry_diagnostics.json').write_text(json.dumps({'baseline':baseline,'changes':diagnostics,'scope':'counts of enumerated and surviving bboxes, not a direct proof of individual newly active faces'},indent=2)+'\n')
 u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'points':len(tasks),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(results)-len(errors),'errors':len(errors),'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(args.out/'results.json'),'top':top};(args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'done={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
 if errors:raise RuntimeError('candidate failure')
if __name__=='__main__':main()
