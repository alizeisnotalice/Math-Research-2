#!/usr/bin/env python3
"""Scan every corresponding-coordinate supercell group and bracket the winner."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
from scan_supercells import canon,digest,make,worker

def run_batch(base,inputpath,out,workers,index,specs,cache,probehash):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True);unique={}
    for spec in specs:
        if spec['kind']=='supercell_baseline':d=make(base,spec['axis'])
        else:d=make(base,spec['axis'],spec['members'],1,F(spec['shift']))
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
    p=argparse.ArgumentParser();p.add_argument('--input',required=True,type=Path);p.add_argument('--probe',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--workers',type=int,default=12);args=p.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    base=json.loads(args.input.read_text());pm=json.loads((args.probe/'manifest.json').read_text());pr=json.loads((args.probe/'results.json').read_text());calcsha=sha(Path(__file__).with_name('periodic_halo.py'))
    if pm['calculator_sha256']!=calcsha or pm['bits']!=40 or pm['initial_tail_radius']!=5 or pm['source_input_sha256']!=sha(args.input):raise ValueError('probe incompatible with current calculator/input/bit/tail rules')
    results={r['id']:r for r in pr if 'error' not in r};cache={}
    for t in pm['tasks']:
        r=results.get(t['id'])
        if r and r['input_sha256']==t['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r)
    specs=[];group_counts=[]
    for axis in (0,1):
        gs={}
        for i,p0 in enumerate(base['pts']):gs.setdefault(F(p0[axis]),[]).append(i)
        group_counts.append(len(gs));specs.append({'kind':'supercell_baseline','axis':axis})
        for v,ids in sorted(gs.items()):
            for eps in (F(1,10000),F(1,100000)):
                for sign in (-1,1):specs.append({'kind':'shift_second_copy','axis':axis,'members':ids,'coordinate':str(v),'shift':str(sign*eps)})
    print(f'groups_x={group_counts[0]} groups_y={group_counts[1]} total={len(group_counts) and sum(group_counts)}',flush=True)
    tasks,res,top=run_batch(base,args.input,args.out,args.workers,0,specs,cache,sha(args.probe/'results.json'));byid={t['id']:t for t in tasks};baselines=[r for r in res if any(a['kind']=='supercell_baseline' for a in byid[r['id']]['aliases'])];winner=top[0]
    if F(winner['inner_ratio'])<=max(F(r['outer_ratio']) for r in baselines):print('no strict improvement; stop',flush=True);return
    alias=next((a for a in winner['aliases'] if a['kind']=='shift_second_copy'),None)
    if alias is None:print('winner baseline; stop',flush=True);return
    s=F(alias['shift']);bracket=[{'kind':'shift_second_copy','axis':alias['axis'],'members':alias['members'],'coordinate':alias['coordinate'],'shift':str(s*F(i,4)),'bracket_index':i} for i in range(-16,17)]
    cache2=cache.copy()
    for t in tasks:
        r=next(x for x in res if x['id']==t['id']);cache2[t['input_sha256']]=(t,r)
    run_batch(base,args.input,args.out,args.workers,1,bracket,cache2,sha(args.out/'batch0'/'results.json'))
if __name__=='__main__':main()
