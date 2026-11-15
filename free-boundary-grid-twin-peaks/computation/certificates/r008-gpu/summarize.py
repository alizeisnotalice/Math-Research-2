#!/usr/bin/env python3
"""Recompute weighted diagnostics; no confidence certificate or iid claim."""
from pathlib import Path
from collections import defaultdict
import json, math
base=Path(__file__).resolve().parent
out={}
for name in ['smc','smc-extension','lshape','lshape-batch2','boundary-family']:
    groups=defaultdict(list)
    for line in (base/name/'smc.jsonl').read_text().splitlines():
        x=json.loads(line)
        if x.get('kind')=='replicate':groups[x['shape']].append(x)
    rows=[]
    for shape,rr in groups.items():
        maxz=max(x['logZ'] for x in rr); w=[math.exp(x['logZ']-maxz) for x in rr]; sw=sum(w)
        p=sum(t*x['p'] for t,x in zip(w,rr))/sw
        boundary=7 if name.startswith('lshape') else (rr[0].get('boundary',6))
        rows.append(dict(shape=shape,replicates=len(rr),boundary=boundary,p_Zweighted=p,scaled=boundary*p,Z_ESS=sw*sw/sum(t*t for t in w),top_Z_share=max(w)/sw,min_roots=min(x['min_roots'] for x in rr),metal_gpu_seconds=sum(x['gpu_seconds'] for x in rr),submit_wait_wall_seconds=sum(x['wall_seconds'] for x in rr)))
    out[name]=rows
print(json.dumps(out,indent=2))
