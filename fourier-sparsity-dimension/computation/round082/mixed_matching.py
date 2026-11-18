"""Verify all effective queries of dimension one or two via exact sufficient bounds."""
import hashlib
import json
import multiprocessing as mp
import platform
import random
import socket
import sys
import time

n=12
size=1<<n
leaves=size//2
degree=leaves+n+1
masks=[]

def margin(a,f):
    M=size>>a
    t=M-2*f
    assert t>=0
    # At most 2-a independent identifications after the pure address queries.
    cost=M*(f+2*t)//2-2*(2-a)*M
    return degree**2*cost-leaves**2*(f+t+n-a+1)**2

def pair_batch(first):
    start=time.process_time()
    low=None
    count=0
    for a in range(first,min(first+16,size)):
        for b in range(a+1,size):
            if b<(a^b):
                for x in masks[a]:
                    for y in masks[b]:
                        m=margin(2,(x&y).bit_count())
                        low=m if low is None else min(low,m)
                        count+=1
    return count,low,time.process_time()-start

if __name__=='__main__':
    seed,workers=map(int,sys.argv[1:])
    start=time.perf_counter()
    xs=list(range(size))
    random.Random(seed).shuffle(xs)
    edges=list(zip(xs[::2],xs[1::2]))
    masks.append(None)
    for a in range(1,size):
        h=[0,0]
        for i,(x,y) in enumerate(edges):
            b=(a&x).bit_count()%2
            if b==(a&y).bit_count()%2:
                h[b]|=1<<i
        masks.append(h)
    minimum1=min(margin(1,x.bit_count()) for h in masks[1:] for x in h)
    with mp.get_context('fork').Pool(workers) as pool:
        rows=list(pool.imap_unordered(pair_batch,range(1,size,16)))
    minimum2=min(m for _,m,_ in rows if m is not None)
    nodata=[degree**2*(2*leaves**2-3*leaves*k)-2*leaves**2*(degree-k)**2 for k in (1,2)]
    assert min(*nodata,minimum1,minimum2)>0
    print(json.dumps(dict(python=platform.python_version(),host=socket.gethostname(),n=n,
        workers=workers,seed=seed,edge_sha256=hashlib.sha256(json.dumps(edges).encode()).hexdigest(),
        leaves=leaves,parent_D=degree,rank_one_address_cells=2*(size-1),
        rank_two_address_cells=sum(r[0] for r in rows),pure_address_rank_zero_margins=nodata,
        rank_one_min_margin=minimum1,rank_two_min_margin=minimum2,
        pair_cpu_seconds=sum(r[2] for r in rows),wall_seconds=time.perf_counter()-start),indent=2))
