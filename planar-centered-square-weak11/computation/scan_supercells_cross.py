#!/usr/bin/env python3
"""Scan every corresponding-coordinate supercell group and bracket the winner."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
from scan_supercells import canon,digest,worker

def make(base,axis,members=None,sign=0,eps=F(0),move_axis=None):
    if move_axis is None:move_axis=axis
    d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};P=[F(base['Px']),F(base['Py'])];d['Px' if axis==0 else 'Py']=str(2*P[axis]);pts=[]
    for copy in (0,1):
        for i,p in enumerate(base['pts']):
            q=[F(v) for v in p];q[axis]+=copy*P[axis]
            if copy==1 and members is not None and i in members:q[move_axis]+=sign*eps*P[move_axis]
            for k in (0,1):q[k]%=P[k]*(2 if k==axis else 1)
            pts.append([str(v) for v in q])
    d['pts']=pts;return d

def run_batch(base,inputpath,out,workers,index,specs,cache,probehash):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True);unique={}
    for spec in specs:
        if spec['kind']=='supercell_baseline':d=make(base,spec['axis'])
        else:d=make(base,spec['axis'],spec['members'],1,F(spec['shift']),spec['move_axis'])
        h=digest(d)
        if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
        elif canon(unique[h]['input'])!=canon(d):raise RuntimeError('hash collision')
        unique[h]['aliases'].append(spec)
    tasks=list(unique.values());reuse={}
    for t in tasks:
        if t['input_sha256'] in cache:
            old,record=cache[t['input_sha256']]
            if canon(t['input'])!=canon(old['input']):raise RuntimeError('hash collision')
            reuse[t['id']]={**record,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':probehash}
    manifest={'schema':'supercell-all-groups-v1','batch':index,'input_sha256':sha(inputpath),'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_sha256':sha(__file__),'probe_results_sha256':probehash,'task_list_sha256':digest(tasks),'bits':40,'tail_rule':'start at 5; increment integer until exact bound passes','workers':workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'batch={index} pid={os.getpid()} points={len(tasks)} reused={len(reuse)} new={len(tasks)-len(reuse)}',flush=True)
    completed=path/'completed.jsonl'
    if completed.exists():raise FileExistsError('existing output')
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with completed.open('x') as stream,cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for r in reuse.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
        futs={pool.submit(worker,t):t for t in tasks if t['id'] not in reuse}
        for f in cf.as_completed(futs):
            try:r=f.result();r['reused']=False
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'batch':index,'points':len(tasks),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(results)-len(errors),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top};(path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'batch={index} completed={len(results)-len(errors)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors:raise RuntimeError('candidate evaluation failed')
    return tasks,results,top

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True,type=Path);p.add_argument('--probe',required=True,type=Path);p.add_argument('--same-axis',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--workers',type=int,default=12);args=p.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    base=json.loads(args.input.read_text());pm=json.loads((args.probe/'manifest.json').read_text());pr=json.loads((args.probe/'results.json').read_text());calcsha=sha(Path(__file__).with_name('periodic_halo.py'))
    if pm['calculator_sha256']!=calcsha or pm['bits']!=40 or pm['initial_tail_radius']!=5 or pm['source_input_sha256']!=sha(args.input):raise ValueError('probe incompatible')
    old={r['id']:r for r in pr if 'error' not in r};cache={}
    for t in pm['tasks']:
        r=old.get(t['id'])
        if r and r['input_sha256']==t['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r)
    specs=[]
    for axis in (0,1):
        move_axis=1-axis;groups={}
        for i,p0 in enumerate(base['pts']):groups.setdefault(F(p0[move_axis]),[]).append(i)
        for v,ids in sorted(groups.items()):
            for eps in (F(1,10000),F(1,100000)):
                for sign in (-1,1):specs.append({'kind':'shift_second_copy','axis':axis,'move_axis':move_axis,'members':ids,'coordinate':str(v),'shift':str(sign*eps)})
    tasks,res,top=run_batch(base,args.input,args.out,args.workers,0,specs,cache,sha(args.probe/'results.json'))
    same=json.loads((args.same_axis/'batch0'/'summary.json').read_text());cross=top[0];prior=same['top'][0]
    # Fraction lower bounds select a rigorous best available witness; compare to both baselines.
    winner=cross if F(cross['inner_ratio'])>F(prior['inner_ratio']) else prior
    if winner is cross:alias=winner['aliases'][0]
    else:
        alias=winner['aliases'][0].copy();alias['move_axis']=alias['axis']
    same_tasks=json.loads((args.same_axis/'batch0'/'manifest.json').read_text())['tasks'];same_results=json.loads((args.same_axis/'batch0'/'results.json').read_text());lookup={t['id']:t for t in same_tasks};baselines=[r for r in same_results if any(a['kind']=='supercell_baseline' for a in lookup[r['id']]['aliases'])]
    if F(winner['inner_ratio'])<=max(F(r['outer_ratio']) for r in baselines):print('no strict improvement; stop',flush=True);return
    s=F(alias['shift']);bracket=[{'kind':'shift_second_copy','axis':alias['axis'],'move_axis':alias['move_axis'],'members':alias['members'],'coordinate':alias['coordinate'],'shift':str(s*F(i,4)),'bracket_index':i} for i in range(-16,17)]
    # Candidate identities are checked independently; cache records remain auditable.
    cache={}
    for t in same_tasks:
        r=next(x for x in same_results if x['id']==t['id']);cache[t['input_sha256']]=(t,r)
    print(f'combined_winner={winner["id"]} axis={alias["axis"]} move_axis={alias["move_axis"]} shift={alias["shift"]}',flush=True)
    run_batch(base,args.input,args.out,args.workers,1,bracket,cache,sha(args.same_axis/'batch0'/'results.json'))
if __name__=='__main__':main()
