"""Test the squared marginal recursion for every four-bit function and direction."""
import argparse
import json
import multiprocessing as mp
import platform
import time
from functools import cache
from math import lcm

def add(space,v):
    return space | frozenset(x^v for x in space)

SPACES={frozenset({0})}
for dim in range(4):
    SPACES |= {add(h,v) for h in SPACES for v in range(16)}
FLATS=[(sum(1<<e for e in h),len(h).bit_length()-1) for h in SPACES]
SCALE=lcm(*range(1,17))

@cache
def cost(s):
    # Grow maximal densest flats after each contraction. The resulting
    # uniform marginals on each layer are feasible by the matroid base theorem.
    s &= ~1
    old,rank=1,0
    value=0
    while s&~old:
        top,num,den,new_rank=old,0,1,rank
        for flat,r in FLATS:
            if r<=rank or flat&old!=old: continue
            count=(s&(flat^old)).bit_count()
            change=r-rank
            if count*den>num*change or (count*den==num*change and r>new_rank):
                top,num,den,new_rank=flat,count,change,r
        value+=den*den*(SCALE//num)
        old,rank=top,new_rank
    return value

@cache
def geometry(s,v):
    tail=[x for x in range(16) if s>>x&1 and (x&v).bit_count()%2]
    base=[]; coord={0:0}
    def extend(x):
        if x in coord: return
        offset=1<<len(base); base.append(x)
        coord.update({y^x:z^offset for y,z in coord.items()})
    for x in tail: extend(x)
    a=len(base)
    for x in range(16):
        if s>>x&1 and (x&v).bit_count()%2==0: extend(x)
    size=1<<a
    active=[x for x in range(1,16) if s>>x&1 and coord[x]>=size]
    rows=[[(x,coord[x]>>a,(-1)**((coord[x]&u).bit_count()%2)) for x in active]
          for u in range(size)] if active else []
    return cost(sum(1<<x for x in tail)),size,rows

def worker(v):
    start=time.process_time()
    masks=[sum(1<<x for x in range(16) if (x&f).bit_count()%2) for f in range(16)]
    checks=nontrivial=tight=0
    for table in range(65536):
        w=[16-2*(table^mask).bit_count() for mask in masks]
        s=sum(1<<x for x,z in enumerate(w) if z)
        if not any((x&v).bit_count()%2 and w[x] for x in range(16)):continue
        local,num,rows=geometry(s,v)
        total=local*num
        for row in rows:
            vals=[0]*16
            for x,b,sign in row: vals[b]+=sign*w[x]
            total+=cost(sum(1<<i for i,val in enumerate(vals) if val))
        required=cost(s)*num
        assert total>=required,(v,table,total,required)
        checks+=1
        nontrivial+=bool(rows)
        tight+=total==required
    return dict(direction=v,cases=checks,with_data=nontrivial,equality=tight,
                cpu_seconds=time.process_time()-start)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=15)
    args=parser.parse_args()
    then=time.perf_counter()
    with mp.get_context('fork').Pool(args.workers) as pool:
        results=pool.map(worker,range(1,16))
    wall=time.perf_counter()-then
    print(json.dumps(dict(host=platform.node(),python=platform.python_version(),workers=args.workers,
                          tables_per_direction=65536,scale=SCALE,wall_seconds=wall,
                          mean_worker_cores=sum(r['cpu_seconds'] for r in results)/wall,
                          cases=sum(r['cases'] for r in results),results=results),indent=2))
