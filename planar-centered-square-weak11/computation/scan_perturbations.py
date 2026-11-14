#!/usr/bin/env python3
"""Frozen, exhaustive one-coordinate rational perturbation batch."""
import concurrent.futures as cf
import hashlib, json, os, platform, resource, sys, time, traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import calculate, sha

INPUT=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/computation/inputs/legacy_r004_transfer_34_27.json')
OUT=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round002')
EPS=(F(1,1000),F(1,10000),F(1,100000))

def canonical(obj): return json.dumps(obj,sort_keys=True,separators=(',',':')).encode()
def digest(obj): return hashlib.sha256(canonical(obj)).hexdigest()

def tasks(data):
    anchor=max(range(len(data['pts'])),key=lambda i:F(data['pts'][i][2]))
    result=[]
    for eps in EPS:
        for i in range(len(data['pts'])):
            if i != anchor:
                for sign in (-1,1): result.append(dict(kind='mass',i=i,anchor=anchor,sign=sign,eps=str(eps)))
        for i in range(len(data['pts'])):
            for axis in (0,1):
                for sign in (-1,1): result.append(dict(kind='coordinate',i=i,axis=axis,sign=sign,eps=str(eps)))
        for axis in (0,1):
            for sign in (-1,1): result.append(dict(kind='period_scale',axis=axis,sign=sign,eps=str(eps)))
    for i,t in enumerate(result): t['id']=f'{i:04d}'
    return result

def transform(data,t):
    d={k:v for k,v in data.items() if k not in ('reach','source','float_ratio')}
    d['pts']=[p.copy() for p in data['pts']]
    eps=F(t['eps']); sign=t['sign']; kind=t['kind']
    if kind=='mass':
        i,a=t['i'],t['anchor']; mi,ma=F(d['pts'][i][2]),F(d['pts'][a][2]); delta=sign*min(mi,ma)*eps
        d['pts'][i][2]=str(mi+delta); d['pts'][a][2]=str(ma-delta)
    elif kind=='coordinate':
        i,axis=t['i'],t['axis']; period=F(d['Px' if axis==0 else 'Py'])
        d['pts'][i][axis]=str((F(d['pts'][i][axis])+sign*period*eps)%period)
    else:
        axis=t['axis']; key='Px' if axis==0 else 'Py'; factor=1+sign*eps
        d[key]=str(F(d[key])*factor)
        for p in d['pts']: p[axis]=str(F(p[axis])*factor)
    return d

def worker(data,t):
    start=time.monotonic()
    try:
        d=transform(data,t); radius=5
        while True:
            try: r=calculate(d,40,F(radius)); break
            except ValueError as e:
                if str(e)!='tail bound does not prove the requested radius cutoff': raise
                radius+=1
        return {'id':t['id'],'transform':t,'input_sha256':digest(d),'tail_radius':radius,'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'inner_ratio_float':r['inner_ratio_float'],'outer_ratio_float':r['outer_ratio_float'],'elapsed_seconds':time.monotonic()-start}
    except Exception:
        return {'id':t['id'],'transform':t,'elapsed_seconds':time.monotonic()-start,'error':traceback.format_exc()}

def main():
    if not os.path.ismount('/Volumes/data'): raise RuntimeError('data is not mounted')
    data=json.loads(INPUT.read_text()); batch=tasks(data)
    manifest={'schema':'periodic-halo-perturbations-v1','task_count':len(batch),'tasks':batch,'task_list_sha256':digest(batch),'source_input_sha256':sha(INPUT),'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_code_sha256':sha(__file__),'workers':12,'bits':40,'tail_radius_initial':5,'transformed_input_hash_format':'SHA256 UTF-8 json.dumps(sort_keys=True,separators=(comma,colon))','python':sys.version,'platform':platform.platform()}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'PID={os.getpid()} tasks={len(batch)} path={OUT}',flush=True)
    start=time.monotonic(); last=start; results=[]
    with cf.ProcessPoolExecutor(max_workers=12) as pool:
        futs=[pool.submit(worker,data,t) for t in batch]
        for f in cf.as_completed(futs):
            results.append(f.result()); now=time.monotonic()
            if now-last>=60:
                print(f'completed={len(results)}/{len(batch)} seconds={now-start:.1f}',flush=True); last=now
    results.sort(key=lambda r:r['id'])
    (OUT/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    good=[r for r in results if 'error' not in r]; errors=[r for r in results if 'error' in r]
    top=sorted(good,key=lambda r:F(r['inner_ratio']),reverse=True)[:8]
    for r in top:
        d=transform(data,r['transform']); assert digest(d)==r['input_sha256']
        (OUT/f"candidate_{r['id']}.json").write_text(json.dumps(d,indent=2)+'\n')
    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    summary={'task_count':len(batch),'completed':len(good),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':usage.ru_utime,'children_system_cpu_seconds':usage.ru_stime,'children_peak_rss_bytes_macos':usage.ru_maxrss,'task_list_sha256':manifest['task_list_sha256'],'results_sha256':sha(OUT/'results.json'),'top':top}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    return bool(errors)
if __name__=='__main__': sys.exit(main())
