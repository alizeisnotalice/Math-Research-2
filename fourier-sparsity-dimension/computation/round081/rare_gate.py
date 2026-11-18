"""Integer verification of the rare-address gate and the affine factor obstruction."""
import json
import multiprocessing as mp
import platform
import socket
import sys
import time

def walsh(row):
    out = list(row)
    w = 1
    while w < len(out):
        for i in range(0,len(out),2*w):
            for j in range(i,i+w):
                a,b = out[j],out[j+w]
                out[j],out[j+w] = a+b,a-b
        w *= 2
    return out

def rank(labels):
    pivots = {}
    for x in labels:
        while x and x.bit_length() in pivots:
            x ^= pivots[x.bit_length()]
        if x:
            pivots[x.bit_length()] = x
    return len(pivots)

def affine_dimension(labels):
    return rank(x ^ labels[0] for x in labels)

def test_batch(params):
    start, stop, n, k = params
    t0 = time.process_time()
    N, Q = 1<<n, 1<<k
    counts = [0,0]
    for mask in range(start,stop):
        row = [1-2*((mask>>j)&1) for j in range(N)]
        f = walsh(row)
        s = [j for j,w in enumerate(f) if w]
        gated = row+[1]*(N*(Q-1))
        spectrum = walsh(gated)
        expected = set()
        for j in set(s)|{0}:
            for a in range(Q):
                idx = a*N+j
                expected.add(idx)
                v = f[j] - (N if j==0 else 0) + (N*Q if idx==0 else 0)
                assert spectrum[idx]==v
        assert {j for j,v in enumerate(spectrum) if v}==expected
        assert rank(expected)==rank(s)+k
        assert spectrum[0]>0
        exceptions = [i for i,a in enumerate(gated) if a<0]
        assert affine_dimension(exceptions)<=n
        counts[0]+=1
        if n==3:
            # Enumerate an independent three-bit factor. Both sign classes must
            # span the product affinely unless both factors are characters up to sign.
            for j in range(1,255):
                signs = [1-2*((j>>h)&1) for h in range(8)]
                both_char = len(s)==1 and sum(bool(v) for v in walsh(signs))==1
                for sign in (1,-1):
                    pts = [8*i+h for i,x in enumerate(row) for h,z in enumerate(signs)
                           if x*z==sign]
                    assert affine_dimension(pts)==n+3-int(both_char)
                counts[1]+=1
    return *counts,time.process_time()-t0

if __name__=='__main__':
    side,workers = sys.argv[1],int(sys.argv[2])
    cases = [(3,3),(4,2)]
    jobs=[]
    for n,k in cases:
        end=(1<<(1<<n))-1
        for lo in range(1,end,64):
            if ((lo-1)//64 % 3 == 0) == (side=='mini'):
                jobs.append((lo,min(end,lo+64),n,k))
    begin=time.perf_counter()
    with mp.get_context('fork').Pool(workers) as pool:
        totals=list(map(sum,zip(*pool.imap_unordered(test_batch,jobs,chunksize=1))))
    print(json.dumps(dict(host=socket.gethostname(),python=platform.python_version(),
        workers=workers,side=side,gates=totals[0],products=totals[1],
        cpu_seconds=totals[2],wall_seconds=time.perf_counter()-begin),indent=2))
