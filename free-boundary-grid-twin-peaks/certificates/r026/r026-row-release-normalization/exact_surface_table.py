#!/usr/bin/env python3
"""Exact n=3..6 top-row surface table and two-row comparison."""
from __future__ import annotations
from fractions import Fraction
from pathlib import Path
import os
import json

CERT_ROOT=Path(os.environ.get('A022_CERTIFICATES_ROOT',Path(__file__).resolve().parent.parent.parent))
HT={3:156,4:4213440,5:153450311313600,6:20779014929894865426585600}
HR={3:6,4:5600,5:3971406600,6:6637177983794054400}

def ladder_connected_orders(n:int)->int:
    m=2*n; full=(1<<m)-1; adj=[0]*m
    for y in range(2):
        for x in range(n):
            v=y*n+x
            for xx,yy in ((x-1,y),(x+1,y),(x,1-y)):
                if 0<=xx<n and 0<=yy<2:adj[v]|=1<<(yy*n+xx)
    # Number of ordered connected-growth words with exactly a prescribed set.
    dp=[0]*(1<<m)
    for v in range(m):dp[1<<v]=1
    for S in range(1,full+1):
        val=dp[S]
        if not val:continue
        boundary=0
        z=S
        while z:
            bit=z&-z; v=bit.bit_length()-1; boundary|=adj[v];z-=bit
        z=boundary&~S
        while z:
            bit=z&-z;dp[S|bit]+=val;z-=bit
    return dp[full]

rows=[]
for n in range(3,7):
    d=json.loads((CERT_ROOT/f'r021-root-transport/n{n}-roots.json').read_text())
    A1,A2=int(d['A1']),int(d['A2']); H=HT[n]
    p=Fraction((1<<(n-1))*H,A1)
    # Exact top-row first-merge arm family from r024.
    F=H*(1<<(n-3))*((1<<(n-1))-2)
    ratio=Fraction(F,A2)
    ladder=ladder_connected_orders(n)
    HRn=HR[n]
    rows.append({'n':n,'A1':str(A1),'A2':str(A2),'H_T':str(H),'p_n':str(p),'p_n_decimal':float(p),
      'F':str(F),'F_over_A2':str(ratio),'F_over_A2_decimal':float(ratio),'A2_over_A1':str(Fraction(A2,A1)),
      'ladder_onepeak_orders_a_R':str(ladder),
      'H_R':str(HRn),
      'top_mass_over_first_two_rows_mass':str(Fraction((1<<(n-1))*H,ladder*HRn))})
out={'scope':'Exact finite values only. The ladder calculation is an independent subset DP for a(R).','rows':rows}
Path(__file__).with_name('exact-surface-table.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
