#!/usr/bin/env python3
"""Small independent scalar checks for the fixed-range adversarial proof review."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json,os,sys
assert os.path.ismount('/Volumes/data')
records=[];shell_count=0
for n,L,p in [(8,1,2),(12,2,3),(20,5,5)]:
    h=F(1,2**n);T=R=F(2**n);m=16*L*L*4**p;q=F(L)/h;k=int(R/h)+1
    assert q.denominator==1 and m<=T and (F(m)/F(1,2**(4*n))).denominator==1
    assert 16*L*L*4**(p+1)>T  # p is maximal for these three examples
    shells=[]
    for s in range(p+1):
        r=L*2**(p+s);M=2**(p+s);K=4**s
        assert K<=M and K*m==16*r*r
        assert 4*r<T and 2*r<R and 2*F(r)/h<k
        assert 2*F(r)-h<2*r
        per_class=F(M,K)*(4*r*h/m)
        assert per_class==F(1,4*q) and q*per_class==F(1,4)
        shells.append({'s':s,'r':r,'M':M,'K':K,'per_class_beta_lower':str(per_class)})
        shell_count+=1
    z=F(p+1,4);D=L-1
    assert 8*z==2*(p+1)
    assert 16*(D+1)**2>F(2**n,2**(2*(p+1)))
    records.append({'n':n,'L':L,'p':p,'mass':m,'q':str(q),'z_lower':str(z),'shells':shells})
assert 16>2  # infeasible-mass branch n=1,L=1, strongest z=0 RHS
source=Path(__file__).resolve()
result={'status':'PASS','python':sys.version,'script_sha256':sha256(source.read_bytes()).hexdigest(),
 'arithmetic':'Fraction exact','seed':None,'tolerance':None,'parameter_cases':len(records),
 'shell_cases':shell_count,'records':records,
 'scope':'scalar bookkeeping only; no subset-family or source-area enumeration; universal proof analytically reviewed'}
out=source.with_suffix('.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
