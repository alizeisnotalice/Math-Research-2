"""Exact matching-direction profile and a successful nested address flag."""
import json
import multiprocessing as mp
import platform
import random
import socket
import sys
import time
from fractions import Fraction

def work(args):
    n, seed = args
    started=time.process_time()
    V=1<<n
    N=V//2
    xs=list(range(V))
    random.Random(seed).shuffle(xs)
    edges=list(zip(xs[::2],xs[1::2]))
    diff=[0]*V
    for x,y in edges:
        diff[x^y]+=1
    step=1
    while step<V:
        for i in range(0,V,step*2):
            for j in range(i,i+step):
                a,b=diff[j],diff[j+step]
                diff[j],diff[j+step]=a+b,a-b
        step*=2
    counts=[(N+v)//4 for v in diff[1:]]
    assert all((N+v)%4==0 for v in diff[1:])
    for v in (1,3,V-1):
        for s in (0,1):
            assert counts[v-1]==sum((x&v).bit_count()%2==s==(y&v).bit_count()%2 for x,y in edges)
    h=(n+1)//2
    final=[0]*(1<<h)
    for x,y in edges:
        if (x^y)>>(n-h)==0:
            final[x>>(n-h)]+=1
    answer=next((a for a,c in enumerate(final) if c==0),None)
    assert answer is not None
    flag=[]
    r=1
    for k in range(1,h+1):
        M=V>>k
        s=answer>>(h-k)
        f=sum(x>>(n-k)==s==y>>(n-k) for x,y in edges)
        t=M-2*f
        assert t>0
        P2=Fraction(((f+t)*r+n-k+1)**2,(M//2)*f+M*t)
        flag.append([k,s,f,t,str(P2)])
    assert flag[-1][2]==0 and Fraction(flag[-1][-1])>Fraction((N*r+n+1)**2,N*N)
    return dict(n=n,seed=seed,hyperplanes=2*(V-1),min_full=min(counts),max_full=max(counts),
                first_moment_numerator=sum(counts),flag=flag,
                cpu_seconds=time.process_time()-started)

if __name__=='__main__':
    workers,seed=map(int,sys.argv[1:])
    begin=time.perf_counter()
    with mp.get_context('fork').Pool(workers) as pool:
        data=pool.map(work,[(16,seed+j) for j in range(workers)])
    print(json.dumps(dict(host=socket.gethostname(),python=platform.python_version(),workers=workers,
        data=data,wall_seconds=time.perf_counter()-begin),indent=2))
