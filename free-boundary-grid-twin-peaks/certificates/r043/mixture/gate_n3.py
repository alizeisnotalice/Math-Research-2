#!/usr/bin/env python3
"""Direct n=3 completion recursion gate for every profile and connected prefix."""
from fractions import Fraction
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

binary, out = map(Path, sys.argv[1:3])
if out.exists(): raise SystemExit("fresh output directory required")
out.mkdir(parents=True)
n = 3
adj = [sum(1 << u for u in range(n*n) if abs(u % n-v % n)+abs(u//n-v//n)==1) for v in range(n*n)]
def region(h): return sum(1 << (y*n+x) for y,c in enumerate(h) for x in range(int(c)))
def neighbor(S, H):
    out = 0
    while S:
        b=S & -S; S -= b; out |= adj[b.bit_length()-1]
    return out & H
def connected(S):
    if not S:return False
    seen=front=S & -S
    while front:
        b=front & -front; front -= b
        more=adj[b.bit_length()-1] & S & ~seen; seen |= more; front |= more
    return seen == S
def roots(H):
    band=1; a=0
    for y in range(n):
        for x in range(band):
            v=y*n+x; b=1<<v
            if H & b and not(adj[v] & ~H): a |= b
    return a
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

total_prefixes=total_words=0
receipts=[]
for h in map("".join,itertools.product("12",repeat=n)):
    H=region(h)
    @lru_cache(None)
    def direct(S):
        if S==H:return 1
        B=neighbor(S,H)&~S
        return sum(direct(S | (1<<v)) for v in range(n*n) if B>>v&1)
    masks=[S for S in range(1<<(n*n)) if S & ~H == 0 and S and connected(S)]
    for S in masks:
        target=out/f"{h}-{S:x}.json"
        p=subprocess.Popen([str(binary),"3",h,str(target),"--region-only","--prefix-mask",f"{S:x}"])
        rc=p.wait(); got=json.loads(target.read_text()) if rc==0 else None
        want=direct(S)
        ok=got is not None and int(got["prefix_completion"])==want
        receipts.append({"curve":h,"prefix_mask":f"{S:x}","pid":p.pid,"returncode":rc,"waited_reaped":True,"matched":ok,"output_sha256":digest(target)})
        if not ok: raise SystemExit(f"completion mismatch {h} {S:x}")
    A=sum(direct(1<<v) for v in range(n*n) if roots(H)>>v&1)
    # Enumerate each allowed-root legal word; its telescoping exact sampler mass is 1/A.
    if A:
        def walk(S, prob):
            nonlocal_words[0] += 1 if S==H else 0
            if S==H:
                assert prob == Fraction(1,A); return
            B=neighbor(S,H)&~S; den=direct(S)
            assert den == sum(direct(S|(1<<v)) for v in range(n*n) if B>>v&1)
            for v in range(n*n):
                if B>>v&1: walk(S|(1<<v),prob*Fraction(direct(S|(1<<v)),den))
        for v in range(n*n):
            if roots(H)>>v&1:
                nonlocal_words=[0]; walk(1<<v,Fraction(direct(1<<v),A)); total_words += nonlocal_words[0]
    total_prefixes += len(masks)

result={"schema":"exact-mixture-n3-direct-gate-v1","profiles":8,"connected_nonempty_prefixes":total_prefixes,
        "cpp_queries":len(receipts),"direct_transition_and_word_probability_checked":True,
        "enumerated_allowed_root_words":total_words,"binary_sha256":digest(binary),"receipts":receipts}
(out/"audit.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
print(json.dumps({k:v for k,v in result.items() if k!="receipts"},sort_keys=True))
