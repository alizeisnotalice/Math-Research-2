#!/usr/bin/env python3
"""Frozen complete-halo position grids and up to three conditional refinements."""
import concurrent.futures as cf
import hashlib,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import calculate,sha

ROOT=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round002')
OUT=ROOT/'position-grid'
BASE=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/computation/inputs/legacy_r004_transfer_34_27.json')
BEST=ROOT/'candidate_0302.json'
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def transform(base,x,y):
    d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};d['pts']=[p.copy() for p in base['pts']];d['pts'][50][:2]=[x,y];return d

def worker(base,t):
    start=time.monotonic()
    try:
        d=transform(base,t['x'],t['y']);rad=5
        while True:
            try:r=calculate(d,40,F(rad));break
            except ValueError as e:
                if str(e)!='tail bound does not prove the requested radius cutoff':raise
                rad+=1
        return {'id':t['id'],'x':t['x'],'y':t['y'],'input_sha256':digest(d),'tail_radius':rad,'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def batch(base,center,step,index,baseline=None):
    path=OUT/f'batch{index}';path.mkdir(parents=True,exist_ok=True)
    px,py=F(base['Px']),F(base['Py']);coords={}
    def add(x,y,label):
        key=(str(x%px),str(y%py));coords.setdefault(key,[]).append(label)
    for i in range(-16,17):
        for j in range(-16,17):add(center[0]+px*i/step,center[1]+py*j/step,f'local:{i}:{j}')
    if index==0:
        for i in range(32):
            for j in range(32):add(px*i/32,py*j/32,f'global:{i}:{j}')
        add(F(baseline['pts'][50][0]),F(baseline['pts'][50][1]),'baseline')
    tasks=[{'id':f'{n:04d}','x':xy[0],'y':xy[1],'origins':orig} for n,(xy,orig) in enumerate(coords.items())]
    manifest={'schema':'position-grid-v1','batch':index,'points':len(tasks),'center':list(map(str,center)),'local_step_denominator':step,'local_index_range':[-16,16],'global_grid':[32,32] if index==0 else None,'source_candidate_sha256':sha(BEST),'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'hash_format':'SHA256 UTF-8 json.dumps(sort_keys=True,separators=(comma,colon))','bits':40,'workers':12,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'batch={index} points={len(tasks)} pid={os.getpid()} manifest={path}/manifest.json',flush=True)
    u0=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic();last=start;results=[]
    with cf.ProcessPoolExecutor(max_workers=12) as pool:
        futs={pool.submit(worker,base,t):t for t in tasks}
        for f in cf.as_completed(futs):
            try:r=f.result()
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
            results.append(r);now=time.monotonic()
            if now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    errors=[r for r in results if 'error' in r];good=[r for r in results if 'error' not in r]
    top=sorted(good,key=lambda r:F(r['inner_ratio']),reverse=True)[:8]
    for r in top:
        d=transform(base,r['x'],r['y']);assert digest(d)==r['input_sha256'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(d,indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN)
    summary={'batch':index,'points':len(tasks),'completed':len(good),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top}
    (path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
    if errors:raise RuntimeError('batch failed; see results.json')
    lookup={r['id']:r for r in results};winner=top[0];center_id=next(t['id'] for t in tasks if 'local:0:0' in t['origins']);center_result=lookup[center_id]
    winner_origins=next(t['origins'] for t in tasks if t['id']==winner['id'])
    edge=any(s.startswith('local:') and any(abs(int(k))==16 for k in s.split(':')[1:]) for s in winner_origins)
    improved=F(winner['inner_ratio'])>F(center_result['outer_ratio'])
    return winner,improved,edge

def main():
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data is not mounted')
    OUT.mkdir(parents=True,exist_ok=True);base=json.loads(BEST.read_text());baseline=json.loads(BASE.read_text());center=tuple(F(v) for v in base['pts'][50][:2]);step=1600
    for index in range(4):
        winner,improved,edge=batch(base,center,step,index,baseline)
        print(f'batch={index} improved={improved} edge={edge}',flush=True)
        if not (improved or edge):break
        center=(F(winner['x']),F(winner['y']));step*=4
if __name__=='__main__':main()
