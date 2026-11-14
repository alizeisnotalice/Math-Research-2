#!/usr/bin/env python3
"""Three simultaneous y-patterns on four explicitly labelled periodic copies."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from periodic_halo_transport import calculate_transport
from scan_supercell_2x2 import canon,digest

def make(base,coeff):
    A,B,C=map(F,coeff);P=[F(base['Px']),F(base['Py'])];d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};d['Px'],d['Py']=str(2*P[0]),str(2*P[1]);pts=[];mx=F(0);group=set(range(27,50,2))
    for x in (0,1):
      for y in (0,1):
       for i,p in enumerate(base['pts']):
        q=[F(v) for v in p];dy=F(0)
        if y==1 and i in group:dy+=A*F(3,400000)*P[1]
        if (x,y)==(1,1) and i==50:dy+=B*F(1,400000)*P[1]
        if y==0 and i in group:dy+=C*(1 if x==1 else -1)*F(1,100000)*P[1]
        mx=max(mx,abs(dy));q[0]=(q[0]+x*P[0])%(2*P[0]);q[1]=(q[1]+y*P[1]+dy)%(2*P[1]);pts.append([str(v) for v in q])
    if mx>F(1,32):raise ValueError('delta cap exceeded')
    d['pts']=pts;return d,mx

def worker(t,base):
    start=time.monotonic()
    try:
        r=calculate_transport(t['input'],base,[2,2],F(t['delta']),40,F(5))
        return {'id':t['id'],'input_sha256':t['input_sha256'],'delta':t['delta'],'tail_bound':r['tail_bound'],'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates']}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def batch(base,inputpath,out,workers,index,coeff_list,cache,selection):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True);unique={}
    for coeff in coeff_list:
        d,delta=make(base,coeff);h=digest(d)
        if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'delta':str(delta),'coefficients':[]}
        elif canon(unique[h]['input'])!=canon(d) or unique[h]['delta']!=str(delta):raise RuntimeError('identity collision')
        unique[h]['coefficients'].append(list(map(str,coeff)))
    tasks=list(unique.values());reuse={}
    for t in tasks:
        if t['input_sha256'] in cache:
            old,r,sourcehash=cache[t['input_sha256']]
            if canon(old['input'])!=canon(t['input']) or old['delta']!=t['delta']:raise RuntimeError('reuse mismatch')
            reuse[t['id']]={**r,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':sourcehash}
    manifest={'schema':'supercell-joint-3pattern-v1','batch':index,'selection':selection,'input_sha256':sha(inputpath),'original_kernel_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'transport_kernel_sha256':sha(Path(__file__).with_name('periodic_halo_transport.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'bits':40,'tail_radius':5,'delta_cap':'1/32','labels':'canonical (copy_x,copy_y,base_index)','copy_factors':[2,2],'workers':workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'batch={index} pid={os.getpid()} points={len(tasks)} reused={len(reuse)} new={len(tasks)-len(reuse)}',flush=True)
    completed=path/'completed.jsonl'
    if completed.exists():raise FileExistsError('existing results')
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with completed.open('x') as stream,cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for r in reuse.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
        futs={pool.submit(worker,t,base):t for t in tasks if t['id'] not in reuse}
        for f in cf.as_completed(futs):
            try:r=f.result();r['reused']=False
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['coefficients']=t['coefficients'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'batch':index,'points':len(tasks),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(results)-len(errors),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top};(path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'batch={index} completed={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors:raise RuntimeError('candidate evaluation failed')
    return tasks,results,top

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    base=json.loads(args.input.read_text());grid=[(F(a,4),F(b,4),F(c,4)) for a,b,c in product(range(9),range(9),range(-4,5))]
    tasks,res,top=batch(base,args.input,args.out,args.workers,0,grid,{},'fixed 9x9x9, A/B=0..2 and C=-1..1, step 1/4');byid={t['id']:t for t in tasks};lookup={r['id']:r for r in res};bycoeff={tuple(F(x) for x in c):lookup[t['id']] for t in tasks for c in t['coefficients']};baseline=bycoeff[(F(1),F(1),F(0))];winner=top[0]
    if F(winner['inner_ratio'])<=F(baseline['outer_ratio']):print('no strict improvement; stop',flush=True);return
    center=tuple(map(F,winner['coefficients'][0]));scores=[]
    for k in range(3):
        diffs=[]
        for sign in (-1,1):
            q=list(center);q[k]+=sign*F(1,4);r=bycoeff.get(tuple(q))
            if r:diffs.append(abs(F(winner['inner_ratio'])-F(r['inner_ratio'])))
        scores.append(max(diffs,default=F(0)))
    axes=sorted(range(3),key=lambda k:(scores[k],-k),reverse=True)[:2];fine=[]
    for i,j in product(range(-16,17),repeat=2):
        v=list(center);v[axes[0]]+=F(i,64);v[axes[1]]+=F(j,64);fine.append(tuple(v))
    h=sha(args.out/'batch0'/'results.json');cache={t['input_sha256']:(t,lookup[t['id']],h) for t in tasks};selection={'center':list(map(str,center)),'axes':axes,'adjacent_coarse_inner_absolute_differences':list(map(str,scores)),'selection_scope':'largest absolute one-step changes near winning point; search heuristic','step':'1/64','range':[-16,16]}
    batch(base,args.input,args.out,args.workers,1,fine,cache,selection)
if __name__=='__main__':main()
