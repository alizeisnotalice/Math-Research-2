#!/usr/bin/env python3
"""2x2 supercell quadrant perturbations using a transport tail proposal kernel."""
import argparse,concurrent.futures as cf,hashlib,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
from periodic_halo_transport import calculate_transport

def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def make(base,axis=None,members=(),shift=F(0)):
    d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};P=[F(base['Px']),F(base['Py'])];d['Px'],d['Py']=str(2*P[0]),str(2*P[1]);pts=[];maxdist=F(0)
    group=set((27,29,31,33,35,37,39,41,43,45,47,49))
    for a in (0,1):
      for b in (0,1):
        for i,p in enumerate(base['pts']):
            q=[F(v) for v in p];dx=F(0);dy=F(3,400000)*P[1] if b==1 and i in group else F(0)
            if (a,b)==(1,1) and i in members:
                if axis==0:dx+=shift*P[0]
                else:dy+=shift*P[1]
            maxdist=max(maxdist,abs(dx),abs(dy));q[0]=(q[0]+a*P[0]+dx)%(2*P[0]);q[1]=(q[1]+b*P[1]+dy)%(2*P[1]);pts.append([str(v) for v in q])
    if maxdist>F(1,32):raise ValueError('displacement exceeds 1/32')
    d['pts']=pts;return d,maxdist

def worker(t,base):
    start=time.monotonic()
    try:
        r=calculate_transport(t['input'],base,[2,2],F(t['delta']),40,F(5))
        return {'id':t['id'],'input_sha256':t['input_sha256'],'delta':t['delta'],'tail_bound':r['tail_bound'],'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates'],'window_atoms_distinct':r['window_atoms_distinct']}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def batch(base,inputpath,out,workers,index,specs,cache):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True);unique={}
    for spec in specs:
        if spec['kind']=='baseline':d,delta=make(base)
        else:d,delta=make(base,spec['axis'],spec['members'],F(spec['shift']))
        h=digest(d)
        if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'delta':str(delta),'aliases':[]}
        elif canon(unique[h]['input'])!=canon(d) or unique[h]['delta']!=str(delta):raise RuntimeError('input/delta collision')
        unique[h]['aliases'].append(spec)
    tasks=list(unique.values());reused={}
    for t in tasks:
        if t['input_sha256'] in cache:
            old,r,sourcehash=cache[t['input_sha256']]
            if canon(old['input'])!=canon(t['input']) or old['delta']!=t['delta']:raise RuntimeError('reuse mismatch')
            reused[t['id']]={**r,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':sourcehash}
    manifest={'schema':'supercell-2x2-transport-probe-v1','batch':index,'input_sha256':sha(inputpath),'original_kernel_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'transport_kernel_sha256':sha(Path(__file__).with_name('periodic_halo_transport.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'bits':40,'radius':5,'delta_cap':'1/32','labels':'canonical (copy_x,copy_y,base_index) nesting in that order','copy_factors':[2,2],'workers':workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'batch={index} pid={os.getpid()} tasks={len(tasks)} reused={len(reused)} new={len(tasks)-len(reused)}',flush=True)
    completed=path/'completed.jsonl'
    if completed.exists():raise FileExistsError('existing results')
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with completed.open('x') as stream,cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for r in reused.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
        futs={pool.submit(worker,t,base):t for t in tasks if t['id'] not in reused}
        for f in cf.as_completed(futs):
            try:r=f.result();r['reused']=False
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if len(results)==1 or now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f} last_id={r["id"]} last_seconds={r.get("elapsed_seconds")}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'batch':index,'points':len(tasks),'reused':len(reused),'new_calls':len(tasks)-len(reused),'completed':len(results)-len(errors),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top};(path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'batch={index} done={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors:raise RuntimeError('candidate evaluation failed')
    return tasks,results,top

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    base=json.loads(args.input.read_text());specs=[{'kind':'baseline'}]
    for axis in (0,1):
        groups={}
        for i,p in enumerate(base['pts']):groups.setdefault(F(p[axis]),[]).append(i)
        for v,members in sorted(groups.items()):
            for eps in (F(1,10000),F(1,100000)):
                for sign in (-1,1):specs.append({'kind':'quadrant_shift','axis':axis,'members':members,'coordinate':str(v),'shift':str(sign*eps)})
    tasks,results,top=batch(base,args.input,args.out,args.workers,0,specs,{})
    winner=top[0];baseline=next(r for r in results if r['id']=='0000')
    if F(winner['inner_ratio'])<=F(baseline['outer_ratio']):print('no strict improvement; stop',flush=True);return
    a=winner['aliases'][0]
    if a['kind']=='baseline':print('winner baseline; stop',flush=True);return
    s=F(a['shift']);bracket=[{'kind':'quadrant_shift','axis':a['axis'],'members':a['members'],'coordinate':a['coordinate'],'shift':str(s*F(i,4)),'bracket_index':i} for i in range(-16,17)]
    sourcehash=sha(args.out/'batch0'/'results.json');lookup={r['id']:r for r in results};cache={t['input_sha256']:(t,lookup[t['id']],sourcehash) for t in tasks}
    batch(base,args.input,args.out,args.workers,1,bracket,cache)
if __name__=='__main__':main()
