#!/usr/bin/env python3
"""Exact n=7 surface ratios from the frozen remote IE receipt and archived A1/A2."""
from fractions import Fraction
from pathlib import Path
import os
import hashlib,json

HERE=Path(__file__).resolve().parent
CERT_ROOT=Path(os.environ.get('A022_CERTIFICATES_ROOT',HERE.parent.parent))
receipt=HERE.parent/'r026-row-release-independent-ie-n7'/'audit-receipt.json'
roots=CERT_ROOT/'r021-root-n7'/'reference'/'n7-d4-root-moments.json'
r=json.loads(receipt.read_text())
d=json.loads(roots.read_text())
H=int(r['n7']['n7-toprow']['result']['hS']); HR=int(r['n7']['n7-tworows']['result']['hS']); A1=int(d['A1']); A2=int(d['A2']); n=7
def ladder(n):
    m=2*n; full=(1<<m)-1; adj=[0]*m
    for y in range(2):
        for x in range(n):
            v=y*n+x
            for xx,yy in ((x-1,y),(x+1,y),(x,1-y)):
                if 0<=xx<n and 0<=yy<2: adj[v]|=1<<(yy*n+xx)
    dp=[0]*(1<<m)
    for v in range(m): dp[1<<v]=1
    for S in range(1,full+1):
        if not dp[S]: continue
        b=0; z=S
        while z:
            x=z&-z; z-=x; b|=adj[x.bit_length()-1]
        z=b&~S
        while z:
            x=z&-z; z-=x; dp[S|x]+=dp[S]
    return dp[full]
aR=ladder(n); p=Fraction((1<<(n-1))*H,A1); F=H*(1<<(n-3))*((1<<(n-1))-2)
ratio=Fraction(F,A2); two=Fraction((1<<(n-1))*H,aR*HR)
out={'scope':'Exact finite n=7 addendum; no asymptotic claim.',
     'remote_receipt_sha256':hashlib.sha256(receipt.read_bytes()).hexdigest(),'A1A2_source_sha256':hashlib.sha256(roots.read_bytes()).hexdigest(),
     'n':n,'H_T':str(H),'H_R':str(HR),'A1':str(A1),'A2':str(A2),'a_2xn':str(aR),
     'p_n':str(p),'p_n_decimal':float(p),'F':str(F),'F_over_A2':str(ratio),'F_over_A2_decimal':float(ratio),
     'A2_over_A1':str(Fraction(A2,A1)),'top_mass_over_first_two_rows_mass':str(two),'top_mass_over_first_two_rows_mass_decimal':float(two)}
Path(__file__).with_name('n7-surface-addendum.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
