#!/usr/bin/env python3
"""Exact integer count of connected growth orders = one-peak labelings."""
import json
from pathlib import Path
import sys

def count(n):
    N=n*n
    nbr=[]
    for v in range(N):
        x,y=v%n,v//n; m=0
        for w in (v-1 if x else None, v+1 if x+1<n else None, v-n if y else None, v+n if y+1<n else None):
            if w is not None:m|=1<<w
        nbr.append(m)
    dp={1<<v:1 for v in range(N)}
    for size in range(1,N):
        nxt={}
        for mask,ways in dp.items():
            for v in range(N):
                if not (mask>>v)&1 and nbr[v]&mask:
                    mm=mask|(1<<v);nxt[mm]=nxt.get(mm,0)+ways
        dp=nxt
    return dp[(1<<N)-1]

out={str(n):count(n) for n in (2,3)}
Path(sys.argv[1]).write_text(json.dumps({"method":"exact Python integer connected-prefix DP","A1":out},indent=2)+"\n")
print(out)
