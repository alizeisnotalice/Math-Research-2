#!/usr/bin/env python3
"""Exact center-prefix feasibility, floating localized halo, exploratory only."""
import argparse, concurrent.futures as cf, hashlib, importlib.util, json
import math, multiprocessing, os, platform, random, time
from fractions import Fraction as F
from pathlib import Path

def init(path):
    global geom
    spec=importlib.util.spec_from_file_location('geom',path)
    geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)

def one(job):
    count,seed=job;rng=random.Random(seed)
    pts=set()
    mode=seed%3
    while len(pts)<count:
        if mode==0:
            q=(rng.randint(-10,10),rng.randint(-10,10))
        elif mode==1:
            q=(rng.choice((-10,10)),rng.randint(-10,10))
            if rng.randrange(2):q=q[::-1]
        else:
            q=(rng.choice((-1,1))*rng.randint(5,10),rng.choice((-1,1))*rng.randint(5,10))
        if q!=(0,0):pts.add(q)
    pts=sorted(pts);weights=[rng.randint(1,20) for _ in pts]
    groups={}
    for p,m in zip(pts,weights):
        d=max(abs(p[0]),abs(p[1]));groups[d]=groups.get(d,0)+m
    scale=None;cum=0
    for d,m in sorted(groups.items()):
        cum+=m;bound=F(4*d*d,400*cum)
        scale=bound if scale is None else min(scale,bound)
    masses=[scale*m for m in weights]
    # Every complete closed distance prefix obeys mass<=4*d^2/20^2.
    cum=F(0)
    for d in sorted(groups):
        cum+=scale*groups[d];assert cum<=F(4*d*d,400)
    total=sum(masses);assert 0<total<=1
    rows=geom.closures(pts,[0]*count);rects=[]
    for row in rows:
        mass=sum(masses[i] for i in row['ids'])
        if row['span_squared']>=400*mass:continue
        a,b,c,d=[v/20 for v in row['bbox']];rad=math.sqrt(float(mass))/2
        u=max(-.5,b-rad);v=min(.5,a+rad);s=max(-.5,d-rad);t=min(.5,c+rad)
        if u<v and s<t:rects.append((u,v,s,t))
    area=geom.sweep_areas(rects,[])[0];assert area>0
    return {'count':count,'seed':seed,'mode':mode,'positions':[[str(F(x,20)),str(F(y,20))] for x,y in pts],
            'masses':list(map(str,masses)),'mass':str(total),'area_float':area,
            'mass_over_area_float':float(total)/area,'closure_count':len(rows),'pid':os.getpid()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workers',type=int,default=3)
    ap.add_argument('--per-count',type=int,default=200);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data');t=time.perf_counter()
    jobs=[(n,740000+n*1000+i) for n in (2,3,4,6,8,12) for i in range(args.per_count)]
    with cf.ProcessPoolExecutor(args.workers,mp_context=multiprocessing.get_context('spawn'),
            initializer=init,initargs=(str(args.geometry),)) as pool:rows=list(pool.map(one,jobs))
    top=sorted(rows,key=lambda x:x['mass_over_area_float'],reverse=True)[:20]
    out={'schema':'r074-reverse-credit-exploratory','scope':'Center nonhalo exact rational prefixes; all clipped halo areas floating. No negative search result proves the conjectured coefficient4.',
         'python':platform.python_version(),'parent_pid':os.getpid(),'worker_pids':sorted(set(r['pid'] for r in rows)),
         'pool_joined':True,'seed_rule':'740000+1000*count+i','requested_workers':args.workers,
         'geometry_sha256':hashlib.sha256(args.geometry.read_bytes()).hexdigest(),
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'count':len(rows),'top20':top,'wall_seconds':time.perf_counter()-t,'all_results':rows}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'count':len(rows),'worker_pids':out['worker_pids'],'top':top[0],
                     'number_above_4_plus_1e6':sum(r['mass_over_area_float']>4+1e-6 for r in rows),'wall_seconds':out['wall_seconds']}))

if __name__=='__main__':main()
