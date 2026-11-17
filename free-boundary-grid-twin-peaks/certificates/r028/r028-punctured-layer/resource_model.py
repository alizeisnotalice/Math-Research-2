#!/usr/bin/env python3
"""Exact layer census for independent subsets outside closed(v)."""
from collections import defaultdict
from math import ceil, factorial
from pathlib import Path
import json

def census(n, v):
    x,y=v%n,v//n
    forbidden={(x,y)}|{(xx,yy) for xx,yy in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)) if 0<=xx<n and 0<=yy<n}
    valid=[]
    for yy in range(n):
        allowed=sum(1<<xx for xx in range(n) if (xx,yy) not in forbidden)
        valid.append([m for m in range(1<<n) if not(m&(m<<1)) and not(m&~allowed)])
    d={(0,0):1}
    for yy in range(n):
        nd=defaultdict(int)
        for (prev,k),c in d.items():
            for m in valid[yy]:
                if not prev&m: nd[m,k+m.bit_count()]+=c
        d=nd
    h=defaultdict(int)
    for (_,k),c in d.items(): h[k]+=c
    pairs=[]
    for k in range(1,max(h)+1):
        p,c=h[k-1],h[k]; cap=1
        while 5*cap<8*p: cap*=2
        pairs.append((k,p,c,cap,8*(p+c)+12*cap))
    peak=max(pairs,key=lambda z:z[-1])
    return {"n":n,"v":v,"xy":[x,y],"forbidden_closed_v":sorted(yy*n+xx for xx,yy in forbidden),
            "allowed_count":n*n-len(forbidden),"histogram":{str(k):h[k] for k in sorted(h)},"independent_sets":sum(h.values()),
            "max_layer":max(h,key=h.get),"max_layer_count":max(h.values()),
            "peak_transition":{"k":peak[0],"previous":peak[1],"current":peak[2],"hash_slots":peak[3],"payload_bytes":peak[4],"payload_gib":peak[4]/2**30},
            "factorial_bound_bits":factorial(n*n-1).bit_length()}

out={"n7_center":census(7,24),"small_representatives":[census(3,4),census(4,5),census(5,12),census(6,21)]}
Path('resource-model.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
