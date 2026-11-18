"""Enumerate all one- and two-query common-address cells for explicit matchings."""
import hashlib
import json
import multiprocessing as mp
import platform
import random
import socket
import sys
import time

def audit(seed):
    start=time.process_time()
    n=9
    N=1<<n
    points=list(range(N))
    random.Random(seed).shuffle(points)
    edges=list(zip(points[::2],points[1::2]))
    masks=[None]
    for a in range(1,N):
        hits=[0,0]
        for i,(x,y) in enumerate(edges):
            s=(a&x).bit_count()%2
            if s==(a&y).bit_count()%2:
                hits[s]|=1<<i
        masks.append(hits)
    samples=[(v.bit_count(),1) for hits in masks[1:] for v in hits]
    # Every two-dimensional dual subspace has three nonzero rows {a,b,a+b}.
    for a in range(1,N):
        for b in range(a+1,N):
            if b < (a^b):
                for x in masks[a]:
                    for y in masks[b]:
                        samples.append(((x&y).bit_count(),2))
    failures=sum(not (0<f<(N>>k)//2) for f,k in samples)
    record=dict(seed=seed,sha256=hashlib.sha256(json.dumps(edges).encode()).hexdigest(),
                cells=len(samples),bad=failures,min_internal=min(f for f,k in samples),
                min_singletons=min((N>>k)-2*f for f,k in samples))
    return record,time.process_time()-start

if __name__=='__main__':
    first,workers=map(int,sys.argv[1:])
    start=time.perf_counter()
    with mp.get_context('fork').Pool(workers) as pool:
        rows=list(pool.imap_unordered(audit,range(first,first+2*workers)))
    records=sorted((v for v,c in rows),key=lambda r:r['seed'])
    assert any(r['bad']==0 for r in records)
    print(json.dumps(dict(host=socket.gethostname(),python=platform.python_version(),
        n=9,workers=workers,records=records,cpu_seconds=sum(c for v,c in rows),
        wall_seconds=time.perf_counter()-start),indent=2))
