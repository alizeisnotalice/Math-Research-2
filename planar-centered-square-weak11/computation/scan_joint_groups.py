#!/usr/bin/env python3
"""Frozen simultaneous rational group directions around an exact fixed input."""
import argparse,concurrent.futures as cf,hashlib,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import calculate,sha

def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def transform(base,directions,coeffs):
    data={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};pts=[[F(v) for v in p] for p in base['pts']];orig=[p.copy() for p in pts];periods=[F(base['Px']),F(base['Py'])];total=sum(p[2] for p in pts);moves=[]
    for d,coef in zip(directions,coeffs):
        g=d['group'];ids=g['members'];delta=coef*d['sign']*F(d['eps'])
        if d['kind']=='coordinate':moves.append((g['axis'],ids,delta))
        elif d['kind']=='mass':
            anchor=d['anchor'];G=sum(orig[i][2] for i in ids);q=min(G,orig[anchor][2]);amount=delta*q
            for i in ids:pts[i][2]+=amount*orig[i][2]/G
            pts[anchor][2]-=amount
        else:raise ValueError('unknown kind')
    if any(p[2]<=0 for p in pts):raise ValueError('nonpositive mass')
    if sum(p[2] for p in pts)!=total:raise ValueError('mass not preserved')
    for axis,ids,delta in moves:
        for i in ids:pts[i][axis]+=periods[axis]*delta
    for p in pts:
        for axis in (0,1):p[axis]%=periods[axis]
    data['pts']=[[str(v) for v in p] for p in pts];return data

def worker(t):
    start=time.monotonic()
    try:
        R=5
        while True:
            try:r=calculate(t['input'],40,F(R));break
            except ValueError as e:
                if str(e)!='tail bound does not prove the requested radius cutoff':raise
                R+=1
        return {'id':t['id'],'input_sha256':t['input_sha256'],'tail_radius':R,'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates']}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def batch(base,inputpath,dirs,out,workers,index,center,step):
    path=out/f'batch{index}';path.mkdir(parents=True,exist_ok=True);unique={};rejected=[]
    for i in range(-16,17):
      for j in range(-16,17):
        coeffs=[center[0]+i*step,center[1]+j*step];alias={'i':i,'j':j,'coefficients':list(map(str,coeffs))}
        try:
            d=transform(base,dirs,coeffs);h=digest(d)
            if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
            elif canonical(unique[h]['input'])!=canonical(d):raise RuntimeError('hash collision')
            unique[h]['aliases'].append(alias)
        except ValueError as e:rejected.append({'alias':alias,'reason':str(e)})
    tasks=list(unique.values());manifest={'schema':'simultaneous-group-grid-v1','batch':index,'center':list(map(str,center)),'step':str(step),'index_range':[-16,16],'input_sha256':sha(inputpath),'direction_summaries_sha256':[x['_summary_sha256'] for x in dirs],'directions':dirs,'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'rejected':rejected,'points':len(tasks),'bits':40,'workers':workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'batch={index} points={len(tasks)} rejected={len(rejected)} pid={os.getpid()}',flush=True)
    completed=path/'completed.jsonl'
    if completed.exists():raise FileExistsError('existing results')
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with completed.open('x') as stream,cf.ProcessPoolExecutor(max_workers=workers) as pool:
        futs={pool.submit(worker,t):t for t in tasks}
        for f in cf.as_completed(futs):
            try:r=f.result()
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if now-last>=60:print(f'batch={index} completed={len(results)}/{len(tasks)} seconds={now-start:.1f}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(path/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(path/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'batch':index,'points':len(tasks),'rejected':len(rejected),'completed':len(results)-len(errors),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(path/'results.json'),'top':top};(path/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if errors:raise RuntimeError('candidate evaluation failed')
    centerid=next(t['id'] for t in tasks if any(a['i']==0 and a['j']==0 for a in t['aliases']));baseline=next(r for r in results if r['id']==centerid);winner=top[0];improved=F(winner['inner_ratio'])>F(baseline['outer_ratio']);winner_coeffs=[F(v) for v in winner['aliases'][0]['coefficients']];print(f'batch={index} improved={improved} winner={winner["id"]} lower={winner["inner_ratio"]}',flush=True);return winner_coeffs,improved

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--mass-summary',required=True,type=Path);ap.add_argument('--coordinate-summary',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--workers',type=int,default=12);args=ap.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    dirs=[]
    for p,kind in ((args.mass_summary,'mass'),(args.coordinate_summary,'coordinate')):
        r=json.loads(p.read_text())['top'][0];d=next(a.copy() for a in r['aliases'] if a['kind']==kind);d['_summary_sha256']=sha(p);dirs.append(d)
    base=json.loads(args.input.read_text());args.out.mkdir(parents=True,exist_ok=True);center=[F(0),F(0)];step=F(1,4)
    for index in range(3):
        center,improved=batch(base,args.input,dirs,args.out,args.workers,index,center,step)
        if not improved:break
        step/=4
if __name__=='__main__':main()
