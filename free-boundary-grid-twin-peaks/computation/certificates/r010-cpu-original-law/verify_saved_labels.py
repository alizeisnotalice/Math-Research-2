#!/usr/bin/env python3
"""Independent verifier for saved descending orders from the CPU sampler."""
import json
import sys

p = json.load(open(sys.argv[1]))
n = p["n"]
N = n*n
for item in p["samples"]:
    order = item["descending_order"]
    assert len(order) == N and set(order) == set(range(N))
    label = [0]*N
    for rank, v in enumerate(order):
        label[v] = N-rank
    peaks=[]
    for v in range(N):
        x,y=v%n,v//n
        ngh=[w for w in (v-1 if x else -1,v+1 if x+1<n else -1,v-n if y else -1,v+n if y+1<n else -1) if w>=0]
        if all(label[v] > label[w] for w in ngh): peaks.append(v)
    assert set(peaks) == set(item["reported_peaks"]), (peaks, item)
    a,b=peaks
    assert abs(a%n-b%n)+abs(a//n-b//n) == item["reported_rho"]
print(json.dumps({"status":"passed","n":n,"samples":len(p["samples"])}))
