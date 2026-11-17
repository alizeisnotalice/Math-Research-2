#!/usr/bin/env python3
"""Exact completion counts for every corner monotone n-vertex path, n=3,4."""
from __future__ import annotations
from functools import lru_cache
from itertools import product
from pathlib import Path
import json

def paths(n):
    for word in product('RD',repeat=n-1):
        x=y=0; out=[0]
        for z in word:
            if z=='R':x+=1
            else:y+=1
            out.append(y*n+x)
        yield ''.join(word),tuple(out)

def completion(n,S):
    N=n*n; init=sum(1<<v for v in S); full=(1<<N)-1
    adj=[]
    for v in range(N):
        x,y=v%n,v//n; mask=0
        for xx,yy in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
            if 0<=xx<n and 0<=yy<n:mask|=1<<(yy*n+xx)
        adj.append(mask)
    @lru_cache(None)
    def h(mask):
        if mask==full:return 1
        ans=0; unseen=full^mask
        while unseen:
            bit=unseen&-unseen;unseen-=bit
            v=bit.bit_length()-1
            if adj[v]&mask:ans+=h(mask|bit)
        return ans
    return h(init)

out={'scope':'Exact integer subset DP; checks only n=3,4, not a general proof.','runs':[]}
for n in (3,4):
    rows=[]
    for word,S in paths(n):rows.append({'word':word,'vertices':list(S),'H':completion(n,S)})
    top=next(x['H'] for x in rows if x['word']=='R'*(n-1))
    out['runs'].append({'n':n,'top_H':top,'minimum_H':min(x['H'] for x in rows),
      'all_ge_top':all(x['H']>=top for x in rows),'rows':rows})
Path(__file__).with_name('monotone-prefix-exact.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
