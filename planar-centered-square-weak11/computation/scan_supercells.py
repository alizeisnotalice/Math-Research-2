#!/usr/bin/env python3
"""Full-halo probes which split one group in the second supercell copy."""
import argparse,concurrent.futures as cf,hashlib,json,os,platform,resource,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import calculate,sha

def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def make(base,axis,member_group=None,sign=0,eps=F(0)):
    d={k:v for k,v in base.items() if k not in ('reach','source','float_ratio')};periods=[F(base['Px']),F(base['Py'])];d['Px' if axis==0 else 'Py']=str(2*periods[axis]);pts=[]
    for copy in (0,1):
        for i,p in enumerate(base['pts']):
            q=[F(v) for v in p];q[axis]+=copy*periods[axis]
            if copy==1 and member_group is not None and i in member_group:q[axis]+=sign*eps*periods[axis]
            q[axis]%=2*periods[axis];pts.append([str(v) for v in q])
    d['pts']=pts
    assert len(pts)==2*len(base['pts']) and sum(F(p[2]) for p in pts)==2*sum(F(p[2]) for p in base['pts'])
    return d

def worker(t):
    start=time.monotonic()
    try:
        radius=5
        while True:
            try:r=calculate(t['input'],40,F(radius));break
            except ValueError as e:
                if str(e)!='tail bound does not prove the requested radius cutoff':raise
                radius+=1
        return {'id':t['id'],'input_sha256':t['input_sha256'],'tail_radius':radius,'inner_ratio':r['inner_ratio'],'outer_ratio':r['outer_ratio'],'elapsed_seconds':time.monotonic()-start,'window_x_coordinates':r['window_x_coordinates'],'window_y_coordinates':r['window_y_coordinates'],'window_atoms_distinct':r['window_atoms_distinct']}
    except Exception:return {'id':t['id'],'error':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start}

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--workers',type=int,default=12);args=p.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    if args.workers<1:raise ValueError('workers must be positive')
    base=json.loads(args.input.read_text());args.out.mkdir(parents=True,exist_ok=True);anchor=max(range(len(base['pts'])),key=lambda i:F(base['pts'][i][2]));unique={};selected=[]
    for axis,seeds in ((0,(2,4,anchor)),(1,(27,2,anchor))):
        groups={}
        for seed in seeds:
            value=F(base['pts'][seed][axis]);members=tuple(i for i,p in enumerate(base['pts']) if F(p[axis])==value);groups.setdefault(members,[]).append(seed)
        if len(groups)<3:raise RuntimeError('selected groups overlap; provide at least three groups per axis')
        for members,source in groups.items():selected.append({'axis':axis,'members':members,'seed_indices':source})
        d=make(base,axis);h=digest(d);unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[{'axis':axis,'kind':'supercell_baseline'}]}
        for members,seeds_ in groups.items():
            for eps in (F(1,100000),F(1,10000)):
                for sign in (-1,1):
                    d=make(base,axis,members,sign,eps);h=digest(d);alias={'axis':axis,'kind':'shift_second_copy','members':members,'seed_indices':seeds_,'eps':str(eps),'sign':sign}
                    if h not in unique:unique[h]={'id':f'{len(unique):04d}','input_sha256':h,'input':d,'aliases':[]}
                    elif canon(unique[h]['input'])!=canon(d):raise RuntimeError('hash collision')
                    unique[h]['aliases'].append(alias)
    tasks=list(unique.values());manifest={'schema':'supercell-probe-v1','source_input_sha256':sha(args.input),'source_atom_count':len(base['pts']),'anchor_index':anchor,'selected_groups':selected,'calculator_sha256':sha(Path(__file__).with_name('periodic_halo.py')),'scan_sha256':sha(__file__),'task_list_sha256':digest(tasks),'bits':40,'initial_tail_radius':5,'workers':args.workers,'python':sys.version,'platform':platform.platform(),'tasks':tasks}
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'pid={os.getpid()} tasks={len(tasks)} groups={selected}',flush=True)
    completed=args.out/'completed.jsonl'
    if completed.exists():raise FileExistsError('existing completed.jsonl')
    start=time.monotonic();last=start;u0=resource.getrusage(resource.RUSAGE_CHILDREN);results=[]
    with completed.open('x') as stream,cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs={pool.submit(worker,t):t for t in tasks}
        for f in cf.as_completed(futs):
            try:r=f.result()
            except Exception:r={'id':futs[f]['id'],'error':traceback.format_exc()}
            results.append(r);stream.write(json.dumps(r,separators=(',',':'))+'\n');stream.flush();now=time.monotonic()
            if len(results)==1 or now-last>=60:print(f'completed={len(results)}/{len(tasks)} seconds={now-start:.1f} last_id={r["id"]} last_seconds={r.get("elapsed_seconds")}',flush=True);last=now
    results.sort(key=lambda r:r['id']);(args.out/'results.json').write_text(json.dumps(results,indent=2)+'\n');errors=[r for r in results if 'error' in r];top=sorted((r for r in results if 'error' not in r),key=lambda r:F(r['inner_ratio']),reverse=True)[:8];byid={t['id']:t for t in tasks}
    for r in top:
        t=byid[r['id']];r['aliases']=t['aliases'];(args.out/f"candidate_{r['id']}.json").write_text(json.dumps(t['input'],indent=2)+'\n')
    u=resource.getrusage(resource.RUSAGE_CHILDREN);summary={'points':len(tasks),'completed':len(results)-len(errors),'errors':len(errors),'error_ids':[r['id'] for r in errors],'wall_seconds':time.monotonic()-start,'children_user_cpu_seconds':u.ru_utime-u0.ru_utime,'children_system_cpu_seconds':u.ru_stime-u0.ru_stime,'children_peak_rss_bytes_macos_cumulative':u.ru_maxrss,'manifest_task_list_sha256':digest(tasks),'results_sha256':sha(args.out/'results.json'),'top':top};(args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
    if errors:raise RuntimeError('candidate evaluation failed')
if __name__=='__main__':main()
