#!/usr/bin/env python3
"""Independent forward-growth check of selected n=5 natural-cut records."""
import json
from math import comb
from pathlib import Path

def adjacencies(n):
    a=[0]*(n*n)
    for y in range(n):
        for x in range(n):
            v=y*n+x
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                X,Y=x+dx,y+dy
                if 0<=X<n and 0<=Y<n:a[v]|=1<<(Y*n+X)
    return a

def reflect(mask,n):
    z=0
    for y in range(n):
        for x in range(n):
            if mask>>(y*n+x)&1:z|=1<<(y*n+n-1-x)
    return z

def forward_root(region,root,adj):
    # Unlike the C++ deletion recurrence this builds legal high-to-low orders.
    level={1<<root:1}
    while True:
        after={}; advanced=False
        for seen,count in level.items():
            if seen==region:
                after[seen]=after.get(seen,0)+count;continue
            bd=0;x=seen
            while x:
                bit=x&-x;x-=bit;bd|=adj[bit.bit_length()-1]
            x=bd&region&~seen
            while x:
                bit=x&-x;x-=bit;after[seen|bit]=after.get(seen|bit,0)+count;advanced=True
        level=after
        if not advanced:return level.get(region,0)

def safe_count(region,n,adj):
    return sum(forward_root(region,y*n,adj) for y in range(n)
               if region>>(y*n)&1 and not (adj[y*n]&~region))

def main():
    n=5;N=25;full=(1<<N)-1;adj=adjacencies(n)
    # Full records are intentionally not needed: select spread-out masks from
    # the C++ ledger after its aggregate result was frozen.
    records=json.loads(Path('allsafe-n5-records.json').read_text())['records']
    positions=[i*(len(records)-1)//9 for i in range(10)]
    checked=[]
    for i in positions:
        r=records[i];A=r['A_mask'];B=full^A
        a=safe_count(A,n,adj);b=safe_count(reflect(B,n),n,adj)
        c=comb(N,A.bit_count())*a*b
        assert (a,b,c)==(int(r['a_left']),int(r['a_right']),int(r['raw_count']))
        checked.append({'index':i,'A_mask':A,'a_left':a,'a_right':b,'raw_count':str(c)})
    out={'status':'passed','method':'independent forward-growth','checked_records':checked}
    Path('verify-selected-forward.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
