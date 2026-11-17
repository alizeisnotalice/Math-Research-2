#!/usr/bin/env python3
"""Independent integer subset DP; does not use peak inclusion-exclusion.

Descending labels reveal vertices of H=grid-v. A birth outside N_G(v)
counts as a surviving peak. A2v has exactly one such birth; Bv also
requires the first revealed vertex to be outside N_G(v).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import time


def count(n, hole):
    verts = [v for v in range(n*n) if v != hole]
    def adjacent(u, v):
        return abs(u//n-v//n)+abs(u % n-v % n) == 1
    neighbors = [sum(1 << j for j, w in enumerate(verts) if adjacent(v,w))
                 for v in verts]
    external = [not adjacent(v,hole) for v in verts]
    size = 1 << len(verts)
    zero, one, top = [0]*size, [0]*size, [0]*size
    zero[0] = 1
    for s in range(size-1):
        rem = (size-1)^s
        while rem:
            bit = rem & -rem
            rem -= bit
            j = bit.bit_length()-1
            t = s | bit
            birth = external[j] and not (neighbors[j] & s)
            if birth:
                one[t] += zero[s]
                if s == 0:
                    top[t] += 1
            else:
                zero[t] += zero[s]
                one[t] += one[s]
                top[t] += top[s]
    return one[-1], top[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reference-dir', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    started = time.monotonic()
    runs, refs = [], []
    for n in (2, 3, 4):
        reference = None
        if n >= 3:
            p = args.reference_dir/f'punctured-table-n{n}.json'
            raw = p.read_bytes()
            reference = json.loads(raw)
            refs.append({'file':p.name,'sha256':hashlib.sha256(raw).hexdigest()})
        # One representative for every D4 orbit, chosen geometrically.
        representatives = sorted({min((x,y),(y,x),(n-1-x,y),(x,n-1-y),
            (n-1-x,n-1-y),(n-1-y,x),(y,n-1-x),(n-1-y,n-1-x))
            for x in range(n) for y in range(n)})
        for x,y in representatives:
            v = x*n+y
            a,b = count(n,v)
            if reference is not None:
                row = reference['per_vertex'][v]
                assert row['v'] == v
                assert (a,b) == (row['A2_v'],row['B_v'])
            else:
                assert (a,b) == (2,2)
            runs.append({'n':n,'v':v,'A2_v':str(a),'B_v':str(b)})
    result = {'status':'passed','method':'descending-label subset birth DP',
        'scope':'all D4 hole classes n=2..4; no asymptotic claim',
        'runs':runs,'reference_files':refs,'python':platform.python_version(),
        'pid':os.getpid(),'wall_seconds':time.monotonic()-started,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
