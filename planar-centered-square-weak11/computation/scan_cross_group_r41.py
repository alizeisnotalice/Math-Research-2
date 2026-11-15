#!/usr/bin/env python3
"""Combine distinct strictly improving R26 coordinate directions; exact labelled transport."""
import concurrent.futures as cf, itertools, json, os, platform, resource, sys, time
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
from periodic_halo_dominated import calculate_dominated
from scan_supercell_2x2 import canon,digest
from scan_supercell_mass import target_delta
from scan_quadrant_groups import moved
P=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11'); src=P/'round026/quadrant-groups'; out=P/'round041/cross-group-probe'
basepath=P/'round008/joint-groups/batch2/candidate_0611.json'; targetpath=P/'round035/axis-wall-refine/candidate_0134.json'
def worker(t,base):
 s=time.monotonic()
 try:
  r=calculate_dominated(t['input'],base,[2,2],F(t['delta']),40,F(8),require_mass_conservation=True)
  return {'id':t['id'],'input_sha256':t['input_sha256'],'delta':t['delta'],'kappa':r['transport_matching']['kappa'],'tail_radius':8,'tail_bound':r['tail_bound'],'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'box_counts':r['counts'],'elapsed_seconds':time.monotonic()-s}
 except Exception:
  import traceback; return {'id':t['id'],'error':traceback.format_exc()}
def main():
 if not os.path.ismount('/Volumes/data'): raise RuntimeError('data not mounted')
 base=json.loads(basepath.read_text()); target=json.loads(targetpath.read_text()); m=json.loads((src/'manifest.json').read_text()); results=json.loads((src/'results.json').read_text()); byid={t['id']:t for t in m['tasks']}; assert sha(basepath)==m['base_sha256']
 baseline=F(results[0]['outer_ratio']); improved=sorted((r for r in results if F(r['inner_ratio'])>baseline),key=lambda r:F(r['inner_ratio']),reverse=True); directions=[]; seen=set()
 for r in improved:
  a=byid[r['id']]['aliases'][0]; g=a['group']; key=(a['quadrant'],g['axis'],tuple(g['members']),F(a['shift'])>0)
  if key not in seen: seen.add(key); directions.append({'source_id':r['id'],'alias':a})
 dirs=[{'quadrants':[q], 'members':g} for g in ([26,27],[32,33]) for q in range(4)]; unique={}
 for direction in dirs:
  for j in range(-12,13):
   d=target
   for q in direction['quadrants']: d=moved(base,d,q,0,direction['members'],F(j,1200000))
   h=digest(d);delta=str(target_delta(base,d))
   if h not in unique: unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'delta':delta,'combinations':[]}
   elif canon(unique[h]['input'])!=canon(d) or unique[h]['delta']!=delta: raise ValueError('collision')
   unique[h]['combinations'].append({'quadrants':direction['quadrants'],'members':direction['members'],'step':str(F(j,1200000))})
 # Put baseline first for reuse regardless of iteration order.
 targethash=digest(target); baseentry=unique.pop(targethash); unique={targethash:baseentry,**unique}
 for i,t in enumerate(unique.values()): t['id']=f'{i:04d}'
 tasks=list(unique.values()); out.mkdir(parents=True,exist_ok=False)
 manifest={'schema':'r41-cross-group-probe-v1','base_sha256':sha(basepath),'target_sha256':sha(targetpath),'source_manifest_sha256':sha(src/'manifest.json'),'source_results_sha256':sha(src/'results.json'),'kernel_sha256':sha(Path(__file__).with_name('periodic_halo_dominated.py')),'script_sha256':sha(__file__),'strict_improving_count':len(improved),'distinct_direction_count':len(directions),'selection':'single quadrant x groups [26,27] and [32,33]; j/1200000 for -12<=j<=12; sampled, no exact wall claim','directions':dirs,'bits':40,'radius':8,'workers':12,'source_r35_summary_sha256':sha(P/'round035/axis-wall-refine/summary.json'),'python':sys.version,'platform':platform.platform(),'task_list_sha256':digest(tasks),'tasks':tasks}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('tasks',len(tasks),'all directions',len(directions),'pid',os.getpid(),flush=True)
 start=time.monotonic(); u0=resource.getrusage(resource.RUSAGE_CHILDREN); rows=[];last=start
 prior=json.loads((P/'round035/axis-wall-refine/summary.json').read_text())['top'][0]; assert prior['input_sha256']==tasks[0]['input_sha256'] and prior['delta']==tasks[0]['delta']; reused={**prior,'id':tasks[0]['id'],'reused':True,'source_results_sha256':sha(P/'round035/axis-wall-refine/results.json')}
 with (out/'completed.jsonl').open('x') as stream,cf.ProcessPoolExecutor(max_workers=12) as pool:
  rows.append(reused);stream.write(json.dumps(reused,separators=(',',':'))+'\n');stream.flush()
  for f in cf.as_completed([pool.submit(worker,t,base) for t in tasks[1:]]):
   r=f.result();rows.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
   if time.monotonic()-last>60:print('completed',len(rows),'/',len(tasks),flush=True);last=time.monotonic()
 rows.sort(key=lambda x:x['id']); (out/'results.json').write_text(json.dumps(rows,indent=2)+'\n'); errs=[r for r in rows if 'error' in r]; top=sorted((r for r in rows if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8]
 for r in top:
  t=tasks[int(r['id'])];r['combinations']=t['combinations'];(out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
 u=resource.getrusage(resource.RUSAGE_CHILDREN); summary={'points':len(tasks),'new_calls':len(tasks)-1,'reused':1,'completed':len(rows)-len(errs),'errors':len(errs),'wall_seconds':time.monotonic()-start,'child_user_cpu_seconds':u.ru_utime-u0.ru_utime,'child_system_cpu_seconds':u.ru_stime-u0.ru_stime,'results_sha256':sha(out/'results.json'),'top':top}; (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('done',len(rows),'errors',len(errs),flush=True)
 if errs: raise RuntimeError('errors')
if __name__=='__main__':main()
