#!/usr/bin/env python3
"""Exact scalar and combinatorial checks, not enumeration of enormous tests."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import hashlib,json,os,sys,platform

def ceil_log2_int(x): return (x-1).bit_length()
def main():
    assert os.path.ismount('/Volumes/data')
    layers=[]; shells=0
    for L in (1,2,3):
      for p in range(5):
        m=16*L*L*4**p;n=ceil_log2_int(m);h=F(1,2**n);R=T=F(2**n);q=F(L)/h
        assert q.denominator==1 and (F(m)*2**(4*n)).denominator==1 and m<=2**n
        records=[]
        for s in range(p+1):
            M=2**(p+s);K=4**s;r=F(L*2**(p+s));t=2*r
            assert K<=M and K*m==16*r*r and 2*t<T and t<R
            assert 2*r/h<R/h+1
            participation=F(comb(M-1,K-1),comb(M,K))
            assert participation==F(K,M)
            a_lower=4*r*h
            per_class=F(M,K)*a_lower/m
            assert per_class==F(1,4*q) and q*per_class==F(1,4)
            # Explicit farthest source is 2r-h; test center h/2 and y=t/2.
            d=max(abs(2*r-h-h/2),t/2)
            witness=(d+t)/2
            assert d<witness<t and F(K*m)/(4*witness*witness)>1
            records.append({'s':s,'M':M,'K':K,'radius':str(t),'class_bound':str(per_class),'all_residue_bound':'1/4'})
            shells+=1
        layers.append({'L':L,'p':p,'n':n,'m':m,'q':int(q),'cost_bound':str(F(p+1,4)),'shells':records})
    # Quantitative range implication, squared to avoid irrational arithmetic.
    range_cases=0
    for D in (F(0),F(1,8),F(1),F(3,2),F(7),F(17,3)):
      L=D.numerator//D.denominator+1
      assert D<L<=D+1
      for n in range(1,41):
        if 16*L*L>2**n:
            assert 16*(D+1)**2>2**n
        else:
            p=0
            while 16*L*L*4**(p+1)<=2**n:p+=1
            z=F(p+1,4)
            exponent=n-2*(p+1)
            rhs=F(2**exponent) if exponent>=0 else F(1,2**(-exponent))
            assert 16*(D+1)**2>rhs
            assert 8*z==2*(p+1)
        range_cases+=1
    source=Path(__file__).resolve()
    result={'status':'PASS','python_version':sys.version,'platform':platform.platform(),
      'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'arithmetic':'fractions.Fraction and exact binomial counts','seed':None,'tolerance':None,
      'layers_checked':len(layers),'shells_checked':shells,'quantitative_range_cases':range_cases,'layers':layers,
      'scope':'Exact selected scalar/binomial certificates; full all-input theorem uses analytic proof. No massive subset or full graph enumeration.'}
    out=source.parent/'results/fixed-range-potential/author-certificate.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','layers_checked','shells_checked','quantitative_range_cases')},indent=2));print(out)
if __name__=='__main__': main()
