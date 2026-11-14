#!/usr/bin/env python3
"""Frozen subsets of distinct strictly improving exact mass increments."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from scan_supercell_mass import canon,digest,target_delta,worker

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);ap.add_argument('--target',type=Path,required=True);ap.add_argument('--previous',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());pm=json.loads((args.previous/'batch0/manifest.json').read_text());prior=json.loads((args.previous/'batch0/results.json').read_text());expected=sha(Path(__file__).with_name('periodic_halo_dominated.py'))
    if pm['dominated_kernel_sha256']!=expected or pm['base_sha256']!=sha(args.base) or pm['target_sha256']!=sha(args.target) or pm['initial_radius']!=5 or pm['bits']!=40:raise ValueError('incompatible prior run')
    tasks0=pm['tasks'];byid={t['id']:t for t in tasks0};baseline=next(r for r in prior if r['id']=='0000');orig=[F(p[2]) for p in target['pts']];rank=sorted((r for r in prior if 'error' not in r and F(r['inner_ratio'])>F(baseline['outer_ratio'])),key=lambda r:F(r['inner_ratio']),reverse=True);vectors=[];seen=set()
    for r in rank:
        t=byid[r['id']];vec=tuple(F(p[2])-w for p,w in zip(t['input']['pts'],orig))
        if vec in seen:continue
        if sum(vec)!=0:raise ValueError('mass vector not conserved')
        seen.add(vec);vectors.append({'source_id':r['id'],'source_input_sha256':r['input_sha256'],'aliases':t['aliases'],'vector':list(map(str,vec)),'source_inner_ratio':r['inner_ratio']})
        if len(vectors)==8:break
    vF=[[F(x) for x in v['vector']] for v in vectors];unique={};rejected=[]
    for flags in product((0,1),repeat=len(vectors)):
        masses=[w+sum((v[i] for v,flag in zip(vF,flags) if flag),F(0)) for i,w in enumerate(orig)];alias={'flags':list(flags)}
        if any(w<=0 for w in masses):rejected.append({'alias':alias,'reason':'nonpositive mass'});continue
        if sum(masses)!=sum(orig):raise ValueError('mass not conserved')
        d={k:v for k,v in target.items()};d['pts']=[[p[0],p[1],str(w)] for p,w in zip(target['pts'],masses)];h=digest(d)
        if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
        elif canon(unique[h]['input'])!=canon(d):raise RuntimeError('hash collision')
        unique[h]['aliases'].append(alias)
    tasks=list(unique.values());cache={}
    for batch in ('batch0','batch1'):
        path=args.previous/batch
        if not path.exists():continue
        m=json.loads((path/'manifest.json').read_text());rr={r['id']:r for r in json.loads((path/'results.json').read_text()) if 'error' not in r}
        if m['dominated_kernel_sha256']!=expected or m['base_sha256']!=sha(args.base) or m['target_sha256']!=sha(args.target):raise ValueError('prior cache mismatch')
        for t in m['tasks']:
            r=rr.get(t['id'])
            if r and t['input_sha256']==r['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r,sha(path/'results.json'))
    reuse={}
    for t in tasks:
        if t['input_sha256'] in cache:
            old,r,h=cache[t['input_sha256']]
            if canon(t['input'])!=canon(old['input']):raise RuntimeError('hash collision')
            reuse[t['id']]={**r,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':h}
    delta=target_delta(base,target);args.out.mkdir(parents=True,exist_ok=True);manifest={'schema':'binary-mass-sums-v1','base_sha256':sha(args.base),'target_sha256':sha(args.target),'prior_manifest_sha256':sha(args.previous/'batch0/manifest.json'),'prior_results_sha256':sha(args.previous/'batch0/results.json'),'dominated_kernel_sha256':expected,'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'selected_vectors':vectors,'selection':'top eight distinct exact mass increment vectors with strict inner > baseline outer','rejected':rejected,'bits':40,'initial_tail_radius':5,'delta':str(delta),'workers':args.workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks};(args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'pid={os.getpid()} selected={len(vectors)} points={len(tasks)} reused={len(reuse)} rejected={len(rejected)}',flush=True)
    completed=args.out/'completed.jsonl'
    if completed.exists():raise FileExistsError('existing output')
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with completed.open('x') as stream,cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for r in reuse.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
        futs={pool.submit(worker,t,base,delta):t for t in tasks if t['id'] not in reuse}
        for f in cf.as_completed(futs):
            try:r=f.result();r['reused']=False
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(args.out/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];lookup={t['id']:t for t in tasks}
    for r in top:
        t=lookup[r['id']];r['aliases']=t['aliases'];(args.out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'points':len(tasks),'selected':len(vectors),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(results)-len(errors),'errors':len(errors),'rejected':len(rejected),'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(args.out/'results.json'),'top':top};(args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'done={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors:raise RuntimeError('candidate failure')
if __name__=='__main__':main()
