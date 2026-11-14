#!/usr/bin/env python3
"""Quantize LP proposals, shrink along ray into exact mass bounds, certify halos."""
import argparse,concurrent.futures as cf,json,os,platform,resource,sys,time,traceback
from decimal import Decimal,ROUND_HALF_EVEN
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha
from scan_supercell_mass import canon,digest,target_delta,worker

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--lp-input',type=Path,required=True);ap.add_argument('--lp-output',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    lp=json.loads(args.lp_input.read_text());out=json.loads(args.lp_output.read_text(),parse_float=Decimal);basepath=Path(lp['base_input_path']);targetpath=Path(lp['target_input_path']);base=json.loads(basepath.read_text());target=json.loads(targetpath.read_text())
    if sha(basepath)!=lp['base_input_sha256'] or sha(targetpath)!=lp['target_input_sha256']:raise ValueError('LP input source mismatch')
    cols=[[F(x) for x in col['unit_mass_increment']] for col in lp['columns']];orig=[F(x) for x in lp['target_masses']];hi=[F(x) for x in lp['upper_mass_bound']];lo=[F(x) for x in lp['lower_mass_bound']]
    if any(w<=l or w>h for w,l,h in zip(orig,lo,hi)):raise ValueError('baseline outside bounds')
    unique={};rejected=[]
    for i,p in enumerate(out['proposals']):
      for ray_scale in (F(1,2),F(3,4),F(1),F(5,4),F(3,2),F(2)):
        trust=p['trust_exact'];alias={'proposal_index':i,'gradient_index':p.get('gradient_index'),'gradient_scale':p.get('gradient_scale'),'budget':p.get('budget'),'trust':trust,'ray_scale':str(ray_scale),'exitflag':p.get('exitflag')}
        try:
            if p.get('exitflag',0)<=0:raise ValueError('LP unsuccessful')
            if len(p['u'])!=len(cols):raise ValueError('wrong vector length')
            u=[F(int((Decimal(str(x))*Decimal(10**9)).to_integral_value(rounding=ROUND_HALF_EVEN)),10**9)*F(trust)*ray_scale for x in p['u']]
            diffs=[sum((u[j]*cols[j][i] for j in range(len(cols))),F(0)) for i in range(len(orig))]
            if sum(diffs)!=0:raise ValueError('column sum not zero')
            alpha=F(1)
            for w,l,h,v in zip(orig,lo,hi,diffs):
                if v>0:alpha=min(alpha,(h-w)/v)
                elif v<0:alpha=min(alpha,(w-l)/(-v))
            if alpha<0:raise ValueError('negative feasibility factor')
            if any(w+alpha*v<=l for w,l,v in zip(orig,lo,diffs)):alpha*=F(999999999999,1000000000000)
            vals=[w+alpha*v for w,v in zip(orig,diffs)]
            if any(not l<w<=h for w,l,h in zip(vals,lo,hi)) or sum(vals)!=sum(orig):raise ValueError('exact feasibility failed')
            alias['feasibility_factor']=str(alpha);alias['quantization']='nearest 1e-9 ties-to-even then exact trust';d={k:v for k,v in target.items()};d['pts']=[[p0[0],p0[1],str(w)] for p0,w in zip(target['pts'],vals)];h=digest(d)
            if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
            elif canon(unique[h]['input'])!=canon(d):raise RuntimeError('hash collision')
            unique[h]['aliases'].append(alias)
        except Exception:rejected.append({'alias':alias,'error':traceback.format_exc()})
    tasks=list(unique.values());delta=target_delta(base,target);args.out.mkdir(parents=True,exist_ok=True);manifest={'schema':'lp-mass-quantized-v1','scope':'numerical LP proposes candidates; exact quantization and ray feasibility','lp_input_sha256':sha(args.lp_input),'lp_output_sha256':sha(args.lp_output),'base_sha256':sha(basepath),'target_sha256':sha(targetpath),'dominated_kernel_sha256':sha(Path(__file__).with_name('periodic_halo_dominated.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'input_proposals':len(out['proposals']),'rejected':rejected,'bits':40,'delta':str(delta),'initial_tail_radius':5,'workers':args.workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks};(args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'pid={os.getpid()} points={len(tasks)} aliases={len(out["proposals"])*6} rejected={len(rejected)}',flush=True)
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with (args.out/'completed.jsonl').open('x') as stream,cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs={pool.submit(worker,t,base,delta):t for t in tasks}
        for f in cf.as_completed(futs):
            try:r=f.result()
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(args.out/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(args.out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'points':len(tasks),'completed':len(results)-len(errors),'errors':len(errors),'rejected':len(rejected),'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(args.out/'results.json'),'top':top};(args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(f'done={len(results)} errors={len(errors)} winner={top[0]["id"]}',flush=True)
    if errors or rejected:raise RuntimeError('candidate errors or rejected LPs')
if __name__=='__main__':main()
