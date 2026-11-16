#!/usr/bin/env python3
"""Exact fixed-position caps; floating complete areas for a small D4 probe.

No negative search result is a general bound. Directed certification is a
separate step if an interesting candidate is found. Old R72 code is read only.
"""
import argparse, concurrent.futures as cf, hashlib, importlib.util, json
import math, multiprocessing, os, platform, time
from fractions import Fraction as F
from pathlib import Path

def init(path):
    global geom
    spec=importlib.util.spec_from_file_location('r072_geometry',path)
    geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)

def float_rects(rows,masses,scale):
    ans=[]
    for row in rows:
        mass=sum(masses[i] for i in row['ids'])
        if row['span_squared']>=mass*scale*scale:continue
        a,b,c,d=[float(F(v,scale)) for v in row['bbox']]
        r=math.sqrt(float(mass))/2
        ans.append((b-r,a+r,d-r,c+r))
    return ans

def clipped(rects,box):
    ans=[]
    for a,b,c,d in rects:
        aa=max(a,box[0]);bb=min(b,box[1]);cc=max(c,box[2]);dd=min(d,box[3])
        if aa<bb and cc<dd:ans.append((aa,bb,cc,dd))
    return ans

def run(job):
    n,kind,a_str=job;a=F(a_str);h=F(n+1,n*n);e=F(1,n*n);w=1-e
    inside=[(i*h,j*h) for i in range(n) for j in range(n)]
    if kind=='corners':
        outside=[(x,y) for x in (-a,w+a) for y in (-a,w+a)]
    else:
        outside=[(-a,w/2),(w+a,w/2),(w/2,-a),(w/2,w+a)]
    pts=inside+outside;scale=math.lcm(n*n,*(v.denominator for p in pts for v in p))
    ints=[tuple(int(v*scale) for v in p) for p in pts]
    rows=geom.closures(ints,[0]*len(ints));res=geom.closures(ints[-4:],[0]*4)
    limits=[]
    for row in rows:
        count=sum(i>=n*n for i in row['ids'])
        if count:
            old=sum(i<n*n for i in row['ids'])*e
            limits.append((max(F(1),F(row['span_squared'],scale*scale))-old)/count)
    cap=min(limits);assert cap>=0
    if cap==0:return {'n':n,'kind':kind,'a':a_str,'cap':'0','records':[],'pid':os.getpid()}
    f=lambda z:float(z)
    C=(-.5,f(w+F(1,2)),-.5,f(w+F(1,2)))
    sides=[(-.5,-f(e)),(1.,f(w+F(1,2)))]
    corners=[(x0,x1,y0,y1) for x0,x1 in sides for y0,y1 in sides]
    threshold=F(1,2)+3*e*e
    records=[]
    for ratio in (F(1,5),F(2,5),F(3,5),F(4,5),F(1)):
        m=cap*ratio;masses=[e]*(n*n)+[m]*4
        full=float_rects(rows,masses,scale);residual=float_rects(res,[m]*4,scale)
        delta=geom.sweep_areas(full,residual)[2]
        deltaC=geom.sweep_areas(clipped(full,C),clipped(residual,C))[2]
        corner=sum(geom.sweep_areas(clipped(full,Q),clipped(residual,Q))[2] for Q in corners)
        assert abs(delta-deltaC)<1e-9,(delta,deltaC)
        maxmass=max(sum(masses[i] for i in row['ids']) for row in rows
                    if row['span_squared']<sum(masses[i] for i in row['ids'])*scale*scale)
        assert maxmass==1
        records.append({'ratio':str(ratio),'outside_mass_each':str(m),'Delta_float':delta,
                        'corner_Delta_float':corner,'corner_threshold':str(threshold),
                        'corner_margin_float':float(threshold)-corner,
                        'full_halo_float':geom.sweep_areas(full,[])[0]})
    return {'n':n,'kind':kind,'a':a_str,'cap':str(cap),'records':records,'pid':os.getpid(),
            'positions':[[str(v) for v in p] for p in pts],
            'closure_count':len(rows)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workers',type=int,default=3)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    jobs=[(n,kind,str(a)) for n in (3,5,9) for kind in ('corners','axes')
          for a in (F(1,16),F(1,10),F(3,20),F(1,4),F(3,8),F(1,2),F(3,4),F(1))]
    with cf.ProcessPoolExecutor(args.workers,mp_context=multiprocessing.get_context('spawn'),
            initializer=init,initargs=(str(args.geometry),)) as pool:
        results=list(pool.map(run,jobs))
    flat=[dict(n=x['n'],kind=x['kind'],a=x['a'],**r) for x in results for r in x['records']]
    result={'schema':'r073-root-small-d4-atomic-probe','scope':'Exact finite-position equal-mass admissibility caps; floating complete areas only. No general conclusion or certified extremum.',
        'python':platform.python_version(),'parent_pid':os.getpid(),
        'worker_pids':sorted(set(x['pid'] for x in results)),
        'requested_workers':args.workers,'pool_joined':True,'geometry_count':len(results),
        'area_evaluations':len(flat),'geometry_sha256':hashlib.sha256(args.geometry.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'max_corner':max(flat,key=lambda x:x['corner_Delta_float']),
        'max_Delta':max(flat,key=lambda x:x['Delta_float']),
        'smallest_corner_margin':min(flat,key=lambda x:x['corner_margin_float']),
        'wall_seconds':time.perf_counter()-start,'results':results}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('geometry_count','area_evaluations','worker_pids','max_corner','max_Delta','smallest_corner_margin','wall_seconds')}))

if __name__=='__main__':main()
