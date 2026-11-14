#!/usr/bin/env python3
"""Exactly rationalize numerical LP proposals and certify their complete halos."""
import argparse,concurrent.futures as cf,hashlib,json,os,platform,resource,sys,time,traceback
from decimal import Decimal,ROUND_HALF_EVEN
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import calculate,sha

def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def transform(base,keys,u,scale):
    d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};d['pts']=[p.copy() for p in base['pts']]
    pts=[[F(v) for v in p] for p in d['pts']];original=[p.copy() for p in pts];period=[F(d['Px']),F(d['Py'])];m0=sum(p[2] for p in pts);pdelta=[F(0),F(0)]
    assert len(keys)==len(u)
    for key,val in zip(keys,u):
        z=F(int((Decimal(str(val))*Decimal(10**9)).to_integral_value(rounding=ROUND_HALF_EVEN)),10**9)*scale
        kind,i,axis,anchor=key
        if kind=='coordinate':pts[i][axis]+=period[axis]*z
        elif kind=='mass':
            delta=min(original[i][2],original[anchor][2])*z;pts[i][2]+=delta;pts[anchor][2]-=delta
        elif kind=='period_scale':pdelta[axis]+=z
        else:raise ValueError(f'unknown key {key}')
    assert sum(p[2] for p in pts)==m0
    if any(p[2]<=0 for p in pts):raise ValueError('nonpositive mass')
    for axis,name in enumerate(('Px','Py')):
        factor=1+pdelta[axis]
        if factor<=0:raise ValueError('nonpositive period')
        per=period[axis]*factor;d[name]=str(per)
        for p in pts:p[axis]=(p[axis]*factor)%per
    d['pts']=[[str(v) for v in p] for p in pts]
    return d

def worker(t):
    started=time.monotonic()
    try:
        radius=5
        while True:
            try:r=calculate(t['input'],40,F(radius));break
            except ValueError as e:
                if str(e)!='tail bound does not prove the requested radius cutoff':raise
                radius+=1
        return {'id':t['id'],'input_sha256':t['input_sha256'],'tail_radius':radius,'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-started}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-started}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--lp-input',required=True,type=Path);ap.add_argument('--lp-output',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data is not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    args.out.mkdir(parents=True,exist_ok=True);base=json.loads(args.input.read_text());lp=json.loads(args.lp_input.read_text());out=json.loads(args.lp_output.read_text(),parse_float=Decimal);unique={};rejected=[]
    for n,p in enumerate(out['proposals']):
        for scale in lp['trust_scales']:
            alias={'proposal_index':n,'gradient_row':p['gradient_row'],'group':p['group'],'budget':p['budget'],'scale':scale,'exitflag':p['exitflag']}
            try:
                if p['exitflag']<=0:raise ValueError('unsuccessful LP proposal')
                d=transform(base,lp['keys'],p['u'],F(scale));h=digest(d)
                if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
                elif canonical(unique[h]['input'])!=canonical(d):raise RuntimeError('hash collision')
                unique[h]['aliases'].append(alias)
            except Exception:rejected.append({'alias':alias,'error':traceback.format_exc()})
    tasks=list(unique.values());manifest={'schema':'joint-lp-proposals-v1','scope':'numerical LP proposes candidates only; exact Fraction transformation and halo certification','rounding':'Decimal parsed LP JSON values; nearest 1e-9 ties-to-even; then multiply trust scale','total_aliases':len(out['proposals'])*len(lp['trust_scales']),'unique_candidates':len(tasks),'rejected_aliases':rejected,'input_sha256':sha(args.input),'lp_input_sha256':sha(args.lp_input),'lp_output_sha256':sha(args.lp_output),'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'bits':40,'workers':args.workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'pid={os.getpid()} aliases={manifest["total_aliases"]} unique={len(tasks)} rejected={len(rejected)}',flush=True)
    start=time.monotonic();last=start;results=[];u0=resource.getrusage(resource.RUSAGE_CHILDREN)
    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs={pool.submit(worker,t):t for t in tasks}
        for f in cf.as_completed(futs):
            try:r=f.result()
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
            results.append(r);now=time.monotonic()
            if now-last>=60:print(f'completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(args.out/'results.json').write_text(json.dumps(results,indent=2)+'\n');good=[r for r in results if 'error' not in r];errors=[r for r in results if 'error' in r];top=sorted(good,key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(args.out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'points':len(tasks),'completed':len(good),'errors':len(errors),'rejected_aliases':len(rejected),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(args.out/'results.json'),'top':top}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
    if errors or rejected:sys.exit(1)
if __name__=='__main__':main()
