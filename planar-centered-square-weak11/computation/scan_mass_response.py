#!/usr/bin/env python3
"""Finite responses and joint discrete grid in a frozen mass-vector family."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from scan_supercell_mass import canon,digest,target_delta,worker

def run(base,target,basepath,targetpath,out,workers,index,coeffs,vectors,cache,selection):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True);orig=[F(p[2]) for p in target['pts']];unique={};rejected=[]
    for coef in coeffs:
        masses=[w+sum((c*v[i] for c,v in zip(coef,vectors)),F(0)) for i,w in enumerate(orig)];alias=list(map(str,coef))
        if min(masses)<=0:rejected.append({'coefficients':alias,'reason':'nonpositive mass'});continue
        if sum(masses)!=sum(orig):raise ValueError('mass not conserved')
        d={k:v for k,v in target.items()};d['pts']=[[p[0],p[1],str(w)] for p,w in zip(target['pts'],masses)];h=digest(d)
        if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'coefficients':[]}
        elif canon(unique[h]['input'])!=canon(d):raise RuntimeError('hash collision')
        unique[h]['coefficients'].append(alias)
    tasks=list(unique.values());reuse={}
    for t in tasks:
        if t['input_sha256'] in cache:
            old,r,sourcehash=cache[t['input_sha256']]
            if canon(t['input'])!=canon(old['input']):raise RuntimeError('reuse mismatch')
            reuse[t['id']]={**r,'id':t['id'],'reused':True,'reused_from_id':old['id'],'reused_results_sha256':sourcehash}
    delta=target_delta(base,target);manifest={'schema':'mass-response-joint-v1','batch':index,'selection':selection,'base_sha256':sha(basepath),'target_sha256':sha(targetpath),'dominated_kernel_sha256':sha(Path(__file__).with_name('periodic_halo_dominated.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'bits':40,'initial_radius':5,'delta':str(delta),'workers':workers,'rejected':rejected,'python':sys.version,'platform':platform.platform(),'tasks':tasks};(path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'batch={index} pid={os.getpid()} points={len(tasks)} reused={len(reuse)} rejected={len(rejected)}',flush=True)
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with (path/'completed.jsonl').open('x') as stream,cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for r in reuse.values():results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush()
        futs={pool.submit(worker,t,base,delta):t for t in tasks if t['id'] not in reuse}
        for f in cf.as_completed(futs):
            try:r=f.result();r['reused']=False
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc(),'reused':False}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['coefficients']=t['coefficients'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'batch':index,'points':len(tasks),'reused':len(reuse),'new_calls':len(tasks)-len(reuse),'completed':len(results)-len(errors),'errors':len(errors),'rejected':len(rejected),'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top};(path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'batch={index} completed={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors:raise RuntimeError('candidate failure')
    return tasks,results,top

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);ap.add_argument('--target',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--amplitudes',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());bm=json.loads((args.binary/'manifest.json').read_text());am=json.loads((args.amplitudes/'manifest.json').read_text());kernel=sha(Path(__file__).with_name('periodic_halo_dominated.py'))
    if any(m['base_sha256']!=sha(args.base) or m['target_sha256']!=sha(args.target) or m['dominated_kernel_sha256']!=kernel for m in (bm,am)):raise ValueError('incompatible evidence')
    vectors=[[F(x) for x in v['vector']] for v in bm['selected_vectors']];winner=json.loads((args.amplitudes/'summary.json').read_text())['top'][0];center=tuple(map(F,winner['aliases'][0]['coefficients']));specs=[center]
    for i in range(len(vectors)):
        for step in (F(-1),F(-1,2),F(1,2),F(1)):
            q=list(center);q[i]+=step;specs.append(tuple(q))
    cache={}
    for path,m in ((args.binary,bm),(args.amplitudes,am)):
        rr={r['id']:r for r in json.loads((path/'results.json').read_text()) if 'error' not in r};h=sha(path/'results.json')
        for t in m['tasks']:
            r=rr.get(t['id'])
            if r and r['input_sha256']==t['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r,h)
    tasks,res,top=run(base,target,args.base,args.target,args.out,args.workers,0,specs,vectors,cache,{'center':list(map(str,center)),'steps':['-1','-1/2','1/2','1'],'source_binary_sha256':sha(args.binary/'manifest.json'),'source_amplitudes_sha256':sha(args.amplitudes/'manifest.json')});byid={t['id']:t for t in tasks};rr={r['id']:r for r in res};bycoef={tuple(map(F,c)):rr[t['id']] for t in tasks for c in t['coefficients']};baseline=bycoef[center];dirs=[]
    for i in range(len(vectors)):
        candidates=[]
        for step in (F(-1),F(-1,2),F(1,2),F(1)):
            q=list(center);q[i]+=step;r=bycoef.get(tuple(q))
            if r and F(r['inner_ratio'])>F(baseline['outer_ratio']):candidates.append((F(r['inner_ratio']),step,r['id']))
        if candidates:
            score,step,rid=max(candidates);dirs.append((score,i,step,rid))
    dirs.sort(reverse=True);dirs=dirs[:6]
    if not dirs:print('no strictly improving direction; stop',flush=True);return
    grid=[]
    for z in product((-1,0,1),repeat=len(dirs)):
        q=list(center)
        for s,(_,i,step,_) in zip(z,dirs):q[i]+=s*step
        grid.append(tuple(q))
    h=sha(args.out/'batch0'/'results.json')
    for t in tasks:cache[t['input_sha256']]=(t,rr[t['id']],h)
    selection={'center':list(map(str,center)),'directions':[{'index':i,'step':str(step),'source_id':rid,'source_lower':str(score)} for score,i,step,rid in dirs],'grid':'-1,0,1 multiples of each selected step','criterion':'top six distinct axes with strict improvement over baseline outer'}
    run(base,target,args.base,args.target,args.out,args.workers,1,grid,vectors,cache,selection)
if __name__=='__main__':main()
