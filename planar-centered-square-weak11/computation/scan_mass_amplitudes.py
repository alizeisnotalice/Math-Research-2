#!/usr/bin/env python3
"""Evaluate 1x/2x amplitudes of eight previously frozen mass columns."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from scan_supercell_mass import canon,digest,target_delta,worker

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);ap.add_argument('--target',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());m=json.loads((args.binary/'manifest.json').read_text());expected=sha(Path(__file__).with_name('periodic_halo_dominated.py'))
    if m['base_sha256']!=sha(args.base) or m['target_sha256']!=sha(args.target) or m['dominated_kernel_sha256']!=expected:raise ValueError('input/kernel mismatch')
    vectors=[[F(v) for v in q['vector']] for q in m['selected_vectors']];orig=[F(p[2]) for p in target['pts']];unique={};rejected=[]
    for coeff in product((1,2),repeat=len(vectors)):
        masses=[w+sum((k*v[i] for k,v in zip(coeff,vectors)),F(0)) for i,w in enumerate(orig)];alias={'coefficients':list(coeff)}
        if min(masses)<=0:rejected.append({'alias':alias,'reason':'nonpositive mass'});continue
        if sum(masses)!=sum(orig):raise ValueError('mass not conserved')
        d={k:v for k,v in target.items()};d['pts']=[[p[0],p[1],str(w)] for p,w in zip(target['pts'],masses)];h=digest(d)
        if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
        elif canon(unique[h]['input'])!=canon(d):raise RuntimeError('hash collision')
        unique[h]['aliases'].append(alias)
    tasks=list(unique.values());lookup={r['id']:r for r in json.loads((args.binary/'results.json').read_text()) if 'error' not in r};cache={t['input_sha256']:(t,lookup[t['id']]) for t in m['tasks'] if t['id'] in lookup};reuse={}
    for t in tasks:
        if t['input_sha256'] in cache:
            old,r=cache[t['input_sha256']]
            if canon(old['input'])!=canon(t['input']):raise RuntimeError('cache collision')
            reuse[t['id']]={**r,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':sha(args.binary/'results.json')}
    args.out.mkdir(parents=True,exist_ok=True);manifest={'schema':'mass-amplitudes-v1','base_sha256':sha(args.base),'target_sha256':sha(args.target),'binary_manifest_sha256':sha(args.binary/'manifest.json'),'dominated_kernel_sha256':expected,'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'selected_vectors':m['selected_vectors'],'grid':'each coefficient in {1,2}; all coefficients simultaneous','rejected':rejected,'bits':40,'initial_radius':5,'python':sys.version,'platform':platform.platform(),'workers':args.workers,'tasks':tasks};(args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'pid={os.getpid()} points={len(tasks)} reused={len(reuse)} rejected={len(rejected)}',flush=True)
    delta=target_delta(base,target);start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with (args.out/'completed.jsonl').open('x') as stream,cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for r in reuse.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
        futs={pool.submit(worker,t,base,delta):t for t in tasks if t['id'] not in reuse}
        for f in cf.as_completed(futs):
            try:r=f.result();r['reused']=False
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(args.out/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(args.out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'points':len(tasks),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(results)-len(errors),'errors':len(errors),'rejected':len(rejected),'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(args.out/'results.json'),'top':top};(args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'done={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors:raise RuntimeError('candidate failure')
if __name__=='__main__':main()
