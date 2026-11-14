#!/usr/bin/env python3
"""Frozen complete-halo two-axis period scaling grids."""
import argparse, concurrent.futures as cf, hashlib, json, os, platform, resource, sys, time, traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import calculate, sha

def canonical(x): return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x): return hashlib.sha256(canonical(x)).hexdigest()
def transform(base,a,b):
    d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')}; d['pts']=[p.copy() for p in base['pts']]
    for axis,(key,factor) in enumerate((('Px',1+F(a)),('Py',1+F(b)))):
        if factor<=0:raise ValueError('scale must be positive')
        d[key]=str(F(d[key])*factor)
        for p in d['pts']:p[axis]=str(F(p[axis])*factor)
    return d

def worker(base,t):
    started=time.monotonic()
    try:
        d=transform(base,t['a'],t['b']);radius=5
        while True:
            try:r=calculate(d,40,F(radius));break
            except ValueError as e:
                if str(e)!='tail bound does not prove the requested radius cutoff':raise
                radius+=1
        return {'id':t['id'],'a':t['a'],'b':t['b'],'input_sha256':digest(d),'tail_radius':radius,'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-started}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-started}

def batch(base,input_path,out,workers,center,step,index):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True)
    coords={}
    def add(a,b,label):coords.setdefault((str(a),str(b)),[]).append(label)
    for i in range(-16,17):
        for j in range(-16,17):add(center[0]+F(i,step),center[1]+F(j,step),f'grid:{i}:{j}')
    if index==0:add(F(0),F(1,10000),'r3_directional_best')
    tasks=[{'id':f'{n:04d}','a':ab[0],'b':ab[1],'origins':labels} for n,(ab,labels) in enumerate(coords.items())]
    manifest={'schema':'period-grid-v1','batch':index,'points':len(tasks),'center':list(map(str,center)),'step_denominator':step,'index_range':[-16,16],'input_sha256':sha(input_path),'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'hash_format':'SHA256 UTF-8 json.dumps(sort_keys=True,separators=(comma,colon))','bits':40,'workers':workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'batch={index} points={len(tasks)} pid={os.getpid()} manifest={path}/manifest.json',flush=True)
    u0=resource.getrusage(resource.RUSAGE_CHILDREN);started=time.monotonic();last=started;results=[]
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        futs={pool.submit(worker,base,t):t for t in tasks}
        for f in cf.as_completed(futs):
            try:r=f.result()
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
            results.append(r);now=time.monotonic()
            if now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-started:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    errors=[r for r in results if 'error' in r];good=[r for r in results if 'error' not in r];top=sorted(good,key=lambda r:F(r['inner_ratio']),reverse=True)[:8]
    for r in top:
        d=transform(base,r['a'],r['b']);assert digest(d)==r['input_sha256'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(d,indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'batch':index,'points':len(tasks),'completed':len(good),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-started,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top}
    (path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
    if errors:raise RuntimeError('batch failed; see results.json')
    center_id=next(t['id'] for t in tasks if 'grid:0:0' in t['origins']);center_result=next(r for r in results if r['id']==center_id);winner=top[0]
    return winner,F(winner['inner_ratio'])>F(center_result['outer_ratio'])

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--workers',type=int,default=12);args=p.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data is not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    args.out.mkdir(parents=True,exist_ok=True);base=json.loads(args.input.read_text());center=(F(0),F(0));step=80000
    for index in range(3):
        winner,improved=batch(base,args.input,args.out,args.workers,center,step,index);print(f'batch={index} improved={improved}',flush=True)
        if not improved:break
        center=(F(winner['a']),F(winner['b']));step*=4
if __name__=='__main__':main()
