SOURCE_SHA='d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710'
#!/usr/bin/env python3
"""Deterministic support enumeration; exact singleton-unit certificates, no floats."""
import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time


def supports(degree, max_terms):
    for size in range(2, min(degree + 1, max_terms) + 1):
        for interior in itertools.combinations(range(1, degree), size - 2):
            yield (0, *interior, degree)


def canonical(I, J):
    if (len(I), I) > (len(J), J) and I[-1] == J[-1]:
        return False
    if math.gcd(*I, *J) != 1:
        return False
    pair = (I, J)
    rev = (tuple(I[-1] - i for i in reversed(I)),
           tuple(J[-1] - j for j in reversed(J)))
    if I[-1] == J[-1] and (len(rev[0]), rev[0]) > (len(rev[1]), rev[1]):
        rev = (rev[1], rev[0])
    # Compare using the same total ordering used for swapping.
    key = lambda p: ((len(p[0]), p[0]), (len(p[1]), p[1]))
    return key(pair) <= key(rev)


def classify(I, J):
    diagonals = {}
    for i in I:
        for j in J:
            diagonals.setdefault(i + j, []).append((i, j))
    unitI, unitJ = {0, I[-1]}, {0, J[-1]}
    for edges in diagonals.values():
        if len(edges) == 1:
            i, j = edges[0]
            unitI.add(i)
            unitJ.add(j)
    collisions = False
    for k in sorted(diagonals):
        edges = diagonals[k]
        if len(edges) > 1:
            collisions = True
            for i, j in edges:
                if i in unitI and j in unitJ:
                    return "excluded", [k, i, j]
    return ("unknown" if collisions else "fair"), None


"""No remote filesystem output: complete per-I ledger chunks streamed to controller."""
import base64, io, multiprocessing, platform
from collections import deque

def jobs():
    for total in (45,46):
        for m in range(1,total//2+1):
            n=total-m
            for index,I in enumerate(supports(m,6)):
                yield (m,n,index,I)

def compute(job):
    m,n,index,I=job; counts=dict(raw=0,canonical=0,excluded=0,fair=0,unknown=0)
    buf=io.BytesIO(); digest=hashlib.sha256(); unknown=[]
    with gzip.GzipFile(fileobj=buf,mode='wb',compresslevel=1,mtime=0) as z:
        for J in supports(n,6):
            counts['raw']+=1
            if not canonical(I,J):continue
            state,witness=classify(I,J)
            line=(json.dumps([I,J,state,witness],separators=(',',':'))+'\n').encode()
            z.write(line);digest.update(line);counts['canonical']+=1;counts[state]+=1
            if state=='unknown':unknown.append([I,J])
    data=buf.getvalue()
    return dict(event='chunk',m=m,n=n,index=index,I=I,counts=counts,unknown=unknown,
                raw_sha256=digest.hexdigest(),gzip_sha256=hashlib.sha256(data).hexdigest(),
                gzip_base64=base64.b64encode(data).decode())

def emit(r):print(json.dumps(r,separators=(',',':')),flush=True)

def run(limit=14400,workers=15,smoke=False):
    begin=time.monotonic(); iterator=iter([(2,4,0,(0,2)),(3,5,1,(0,1,3))] if smoke else jobs())
    submitted=finished=0; exhausted=False; pending=[]
    with multiprocessing.get_context('fork').Pool(workers) as pool:
        emit(dict(event='started',pid=os.getpid(),worker_pids=[p.pid for p in pool._pool],
                  workers=workers,python=sys.version,host=platform.node(),source_sha256=SOURCE_SHA,
                  total_min=45,total_max=46,cap=6,seconds_limit=limit,smoke=smoke,
                  storage='RAM only; ledger streamed to local data disk'))
        while True:
            while not exhausted and len(pending)<2*workers and time.monotonic()-begin<limit:
                try:job=next(iterator)
                except StopIteration:exhausted=True;break
                pending.append(pool.apply_async(compute,(job,)));submitted+=1
            ready=[p for p in pending if p.ready()]
            for p in ready:
                emit(p.get());pending.remove(p);finished+=1
            if not pending and (exhausted or time.monotonic()-begin>=limit):break
            if not ready:time.sleep(.01)
        pool.close();pool.join()
    emit(dict(event='complete' if exhausted else 'partial_timeout',workers_joined=True,
              submitted=submitted,finished=finished,seconds=time.monotonic()-begin,
              complete=exhausted,scope='Exact unit/singleton classification only; unresolved cases remain unknown; independent verification still required'))

run(smoke=True,workers=2)
