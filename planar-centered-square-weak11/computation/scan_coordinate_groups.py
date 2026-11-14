#!/usr/bin/env python3
"""Exact whole-coordinate-group translations and mass transfers."""
import argparse,concurrent.futures as cf,hashlib,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import calculate,sha
EPS=('1/100','1/1000','1/10000','1/100000')
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def groups(d):
    answer=[]
    for axis in (0,1):
        m={}
        for i,p in enumerate(d['pts']):m.setdefault(F(p[axis]),[]).append(i)
        answer += [{'axis':axis,'coordinate':str(v),'members':ids} for v,ids in sorted(m.items())]
    return answer

def changed(base,g,kind,sign,eps):
    d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};d['pts']=[p.copy() for p in base['pts']];pts=[[F(v) for v in p] for p in d['pts']];axis=g['axis'];ids=g['members'];e=F(eps)*sign;anchor=None
    if kind=='coordinate':
        P=F(d['Px' if axis==0 else 'Py'])
        for i in ids:pts[i][axis]=(pts[i][axis]+P*e)%P
    else:
        others=[i for i in range(len(pts)) if i not in ids]
        if not others:raise ValueError('group contains all atoms')
        anchor=max(others,key=lambda i:pts[i][2]);G=sum(pts[i][2] for i in ids);q=min(G,pts[anchor][2]);delta=e*q
        for i in ids:pts[i][2]+=delta*pts[i][2]/G
        pts[anchor][2]-=delta
        assert sum(p[2] for p in pts)==sum(F(p[2]) for p in base['pts'])
        if any(p[2]<=0 for p in pts):raise ValueError('nonpositive mass')
    d['pts']=[[str(v) for v in p] for p in pts]
    return d,anchor

def worker(t):
    started=time.monotonic()
    try:
        radius=5
        while True:
            try:r=calculate(t['input'],40,F(radius));break
            except ValueError as e:
                if str(e)!='tail bound does not prove the requested radius cutoff':raise
                radius+=1
        return {'id':t['id'],'input_sha256':t['input_sha256'],'tail_radius':radius,'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-started,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates']}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-started}

def batch(base,source,out,workers,index):
    path=out/f'epoch{index}';path.mkdir(parents=True,exist_ok=True);unique={};skipped=[]
    def add(d,alias):
        h=digest(d)
        if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
        elif canonical(unique[h]['input'])!=canonical(d):raise RuntimeError('hash collision')
        unique[h]['aliases'].append(alias)
    normalized={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};normalized['pts']=[[str(F(v)) for v in p] for p in base['pts']];normalized['Px']=str(F(base['Px']));normalized['Py']=str(F(base['Py']));normalized['lam']=str(F(base['lam']));add(normalized,{'kind':'baseline_and_current_best'})
    gs=groups(normalized)
    for g in gs:
      for kind in ('coordinate','mass'):
       if kind=='mass' and len(g['members'])==len(normalized['pts']):skipped.append({'kind':kind,'group':g,'reason':'no outside anchor'});continue
       for eps in EPS:
        for sign in (-1,1):
            d,anchor=changed(normalized,g,kind,sign,eps);add(d,{'kind':kind,'group':g,'sign':sign,'eps':eps,'anchor':anchor})
    tasks=list(unique.values());manifest={'schema':'coordinate-group-search-v1','epoch':index,'baseline_file_sha256':sha(source),'baseline_canonical_sha256':digest(normalized),'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'eps':EPS,'groups':gs,'skipped':skipped,'unique_candidates':len(tasks),'bits':40,'workers':workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'epoch={index} groups={len(gs)} points={len(tasks)} pid={os.getpid()}',flush=True)
    resultpath=path/'completed.jsonl'
    if resultpath.exists():raise FileExistsError('existing completed.jsonl; preserve existing run')
    start=time.monotonic();last=start;results=[];u0=resource.getrusage(resource.RUSAGE_CHILDREN)
    with resultpath.open('x') as stream,cf.ProcessPoolExecutor(max_workers=workers) as pool:
        futs={pool.submit(worker,t):t for t in tasks}
        for f in cf.as_completed(futs):
            try:r=f.result()
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'epoch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];good=[r for r in results if 'error' not in r];top=sorted(good,key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        task=byid[r['id']];r['aliases']=task['aliases'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(task['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'epoch':index,'points':len(tasks),'completed':len(good),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top}
    (path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
    if errors:raise RuntimeError('candidate evaluation failed')
    baseline=next(r for r in results if r['id']=='0000');winner=top[0];improved=F(winner['inner_ratio'])>F(baseline['outer_ratio']);return path/f"candidate_{winner['id']}.json",improved

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data is not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    args.out.mkdir(parents=True,exist_ok=True);source=args.input
    for index in range(8):
        best,improved=batch(json.loads(source.read_text()),source,args.out,args.workers,index);print(f'epoch={index} improved={improved}',flush=True)
        if not improved:break
        source=best
if __name__=='__main__':main()
