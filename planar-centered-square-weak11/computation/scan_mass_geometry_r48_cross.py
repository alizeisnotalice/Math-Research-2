#!/usr/bin/env python3
"""Cross two previous geometry directions with three exact conserved masses."""
import concurrent.futures as cf,json,os,platform,resource,sys,time
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
from scan_supercell_2x2 import canon,digest
from scan_supercell_mass import target_delta
from scan_quadrant_groups import moved,worker
P=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11');O=P/'round048/mass-geometry';basepath=P/'round008/joint-groups/batch2/candidate_0611.json';targetpath=P/'round043/cross-group-refine/candidate_0195.json'
def main():
 if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
 base=json.loads(basepath.read_text());cache={};bybatch=[]
 for i in range(3):
  p=O/f'batch{i}';m=json.loads((p/'manifest.json').read_text());rs={r['id']:r for r in json.loads((p/'results.json').read_text())};bybatch.append(m)
  for t in m['tasks']:
   r=rs[t['id']];assert digest(t['input'])==t['input_sha256']==r['input_sha256'];cache[t['input_sha256']]=(t,r,sha(p/'results.json'))
 basetarget=json.loads(targetpath.read_text()); chosen=[('baseline',basetarget),('first_grid_winner',json.loads((O/'batch0/candidate_0072.json').read_text())),('final_grid_winner',json.loads((O/'batch2/candidate_0072.json').read_text()))];dirs=[('baseline',[],[],F(0)),('q2_28_minus',[2],[28,29],-F(1,500000)),('q2_28_plus',[2],[28,29],F(1,500000)),('all_26_minus',[0,1,2,3],[26,27],-F(1,500000)),('all_26_plus',[0,1,2,3],[26,27],F(1,500000))];unique={}
 for name,target in chosen:
  for label,qs,ids,shift in dirs:
   d=target
   for q in qs:d=moved(base,d,q,0,ids,shift)
   assert all(F(x[2])>0 for x in d['pts']) and sum(F(x[2]) for x in d['pts'])==sum(F(x[2]) for x in basetarget['pts'])
   h=digest(d);delta=str(target_delta(base,d))
   if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'delta':delta,'aliases':[]}
   elif canon(unique[h]['input'])!=canon(d) or unique[h]['delta']!=delta:raise RuntimeError('hash mismatch')
   unique[h]['aliases'].append({'mass':name,'geometry':label,'shift':str(shift)})
 tasks=list(unique.values());reuse={}
 for t in tasks:
  if t['input_sha256'] in cache:
   old,r,h=cache[t['input_sha256']];assert canon(old['input'])==canon(t['input']);reuse[t['id']]={**r,'id':t['id'],'reused':True,'source_results_sha256':h}
 path=O/'cross';path.mkdir(exist_ok=False);manifest={'schema':'r48-mass-geometry-cross-v1','source_batch_results_sha256':[sha(O/f'batch{i}/results.json') for i in range(3)],'base_sha256':sha(basepath),'target_sha256':sha(targetpath),'kernel_sha256':sha(Path(__file__).with_name('periodic_halo_dominated.py')),'scan_sha256':sha(__file__),'tasks_sha256':digest(tasks),'bits':40,'initial_R':5,'workers':10,'python':sys.version,'platform':platform.platform(),'tasks':tasks};(path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('pid',os.getpid(),'points',len(tasks),'reused',len(reuse),flush=True);start=time.monotonic();u=resource.getrusage(resource.RUSAGE_CHILDREN);rows=[]
 with (path/'completed.jsonl').open('x') as f,cf.ProcessPoolExecutor(max_workers=10) as pool:
  for r in reuse.values():rows.append(r);f.write(json.dumps(r,separators=(',',':'))+'\n');f.flush()
  for future in cf.as_completed([pool.submit(worker,t,base) for t in tasks if t['id'] not in reuse]):
   r=future.result();r['reused']=False;rows.append(r);f.write(json.dumps(r,separators=(',',':'))+'\n');f.flush()
 rows.sort(key=lambda x:x['id']);(path/'results.json').write_text(json.dumps(rows,indent=2)+'\n');errs=[r for r in rows if 'error'in r];top=sorted((r for r in rows if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True);byid={t['id']:t for t in tasks}
 for r in top[:5]:r['aliases']=byid[r['id']]['aliases'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(byid[r['id']]['input'],indent=2)+'\n')
 v=resource.getrusage(resource.RUSAGE_CHILDREN);s={'points':len(tasks),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(rows)-len(errs),'errors':len(errs),'wall_seconds':time.monotonic()-start,'child_user_cpu_seconds':v.ru_utime-u.ru_utime,'results_sha256':sha(path/'results.json'),'top':top[:5]};(path/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print('done',len(rows),'errors',len(errs),'winner',top[0]['id'],flush=True)
 if errs:raise RuntimeError('candidate failure')
if __name__=='__main__':main()
