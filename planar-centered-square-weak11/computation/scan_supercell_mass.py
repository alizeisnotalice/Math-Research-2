#!/usr/bin/env python3
"""Conservative within-quadrant group transfers, dominated-tail proposals."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from periodic_halo_dominated import calculate_dominated
from scan_supercell_2x2 import canon,digest

def apply(base,target,operations):
    d={k:v for k,v in target.items()};pts=[[F(v) for v in p] for p in target['pts']];orig=[p.copy() for p in pts];n=len(base['pts']);total=sum(p[2] for p in pts)
    for op,coef in operations:
        ids=[op['quadrant']*n+i for i in op['members']];anchor=op['quadrant']*n+op['anchor'];G=sum(orig[i][2] for i in ids);q=min(G,orig[anchor][2]);delta=F(coef)*op['sign']*F(op['eps'])*q
        for i in ids:pts[i][2]+=delta*orig[i][2]/G
        pts[anchor][2]-=delta
    if any(p[2]<=0 for p in pts):raise ValueError('nonpositive mass')
    if sum(p[2] for p in pts)!=total:raise ValueError('mass not conserved')
    d['pts']=[[str(v) for v in p] for p in pts];return d

def target_delta(base,target):
    P=[F(base['Px']),F(base['Py'])];n=len(base['pts']);maximum=F(0)
    for a,b,i in product(range(2),range(2),range(n)):
        point=target['pts'][(a*2+b)*n+i];orig=base['pts'][i]
        for axis,k in ((0,a),(1,b)):
            diff=(F(point[axis])-F(orig[axis])-k*P[axis])%(2*P[axis]);maximum=max(maximum,min(diff,2*P[axis]-diff))
    return maximum

def worker(t,base,delta):
    start=time.monotonic();R=5
    try:
        while True:
            try:r=calculate_dominated(t['input'],base,[2,2],delta,40,F(R),require_mass_conservation=True);break
            except ValueError as e:
                if str(e)!='transport tail bound does not prove the requested radius cutoff':raise
                R+=1
        return {'id':t['id'],'input_sha256':t['input_sha256'],'tail_radius':R,'kappa':r['transport_matching']['kappa'],'delta':str(delta),'tail_bound':r['tail_bound'],'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates']}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def batch(base,target,basepath,targetpath,out,workers,index,specs,cache,selection):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True);unique={};reject=[];delta=target_delta(base,target)
    for spec in specs:
        try:
            d=apply(base,target,spec['operations']);h=digest(d)
            if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
            elif canon(unique[h]['input'])!=canon(d):raise RuntimeError('hash collision')
            unique[h]['aliases'].append(spec['alias'])
        except ValueError as e:reject.append({'alias':spec['alias'],'reason':str(e)})
    tasks=list(unique.values());reuse={}
    for t in tasks:
        if t['input_sha256'] in cache:
            old,r,sourcehash=cache[t['input_sha256']]
            if canon(old['input'])!=canon(t['input']):raise RuntimeError('reuse mismatch')
            reuse[t['id']]={**r,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':sourcehash}
    manifest={'schema':'supercell-mass-groups-v1','batch':index,'selection':selection,'base_sha256':sha(basepath),'target_sha256':sha(targetpath),'dominated_kernel_sha256':sha(Path(__file__).with_name('periodic_halo_dominated.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'bits':40,'initial_radius':5,'delta':str(delta),'copy_factors':[2,2],'labels':'canonical (copy_x,copy_y,base_index)','require_mass_conservation':True,'rejected':reject,'workers':workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks};(path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'batch={index} pid={os.getpid()} points={len(tasks)} reused={len(reuse)} rejected={len(reject)}',flush=True)
    resultpath=path/'completed.jsonl'
    if resultpath.exists():raise FileExistsError('existing results')
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with resultpath.open('x') as stream,cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for r in reuse.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
        futs={pool.submit(worker,t,base,delta):t for t in tasks if t['id'] not in reuse}
        for f in cf.as_completed(futs):
            try:r=f.result();r['reused']=False
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'batch':index,'points':len(tasks),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(results)-len(errors),'errors':len(errors),'error_ids':[r['id'] for r in errors],'rejected':len(reject),'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top};(path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'batch={index} completed={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors or reject:raise RuntimeError('candidate failure or rejected input')
    return tasks,results,top

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',required=True,type=Path);ap.add_argument('--target',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());n=len(base['pts']);specs=[{'operations':[],'alias':{'kind':'baseline'}}];skipped=[]
    for quadrant in range(4):
     for axis in (0,1):
      groups={}
      for i,p in enumerate(base['pts']):groups.setdefault(F(p[axis]),[]).append(i)
      for value,ids in sorted(groups.items()):
       others=[i for i in range(n) if i not in ids]
       if not others:skipped.append({'quadrant':quadrant,'axis':axis,'members':ids});continue
       anchor=max(others,key=lambda i:F(target['pts'][quadrant*n+i][2]))
       for eps in ('1/1000','1/10000'):
        for sign in (-1,1):
         op={'quadrant':quadrant,'axis':axis,'members':ids,'coordinate':str(value),'anchor':anchor,'eps':eps,'sign':sign};specs.append({'operations':[(op,'1')],'alias':op})
    tasks,res,top=batch(base,target,args.base,args.target,args.out,args.workers,0,specs,{}, {'skipped_groups':skipped});byid={t['id']:t for t in tasks};lookup={r['id']:r for r in res};baseline=lookup['0000'];rank=sorted((r for r in res if r['id']!='0000' and F(r['inner_ratio'])>F(baseline['outer_ratio'])),key=lambda r:F(r['inner_ratio']),reverse=True);chosen=[];seen=set()
    for r in rank:
        for a in byid[r['id']]['aliases']:
            key=(a['quadrant'],tuple(a['members']),a['anchor'])
            if key not in seen:chosen.append(a);seen.add(key);break
        if len(chosen)==2:break
    if len(chosen)<2:print('fewer than two distinct improving directions; stop',flush=True);return
    fine=[]
    for i,j in product(range(-16,17),repeat=2):
        co=[str(F(i,4)),str(F(j,4))];fine.append({'operations':[(chosen[k],co[k]) for k in range(2)],'alias':{'coefficients':co}})
    h=sha(args.out/'batch0'/'results.json');cache={t['input_sha256']:(t,lookup[t['id']],h) for t in tasks};batch(base,target,args.base,args.target,args.out,args.workers,1,fine,cache,{'directions':chosen,'coefficients':'i/4,j/4 for i,j in [-16,16]','selection_scope':'two highest strict improvements from distinct mass groups'})
if __name__=='__main__':main()
