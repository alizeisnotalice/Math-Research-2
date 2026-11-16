#!/usr/bin/env python3
"""Exact integer/Sturm feasibility audit for the next H_d priority target."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import sympy as S

H=[1,1,1,0,0,0,1,1]
def conv(a,b):
 z=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]+=x*y
 return z
def stretch(a,q):z=[0]*((len(a)-1)*q+1);z[::q]=a;return z
def main():
 x=S.symbols('x');f=sum(c*x**i for i,c in enumerate(H));C=conv(stretch(H,3),stretch(H,4));support=[i for i,c in enumerate(H) if c]
 out={'mother_coefficients_ascending':H,'mother_sha256':hashlib.sha256(json.dumps(H,separators=(',',':')).encode()).hexdigest(),'a':3,'b':4,'degree':len(C)-1,'integer_binary_C':set(C)<={0,1},'C_sha256':hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest(),'support':support,'difference_set_positive':sorted({abs(i-j) for i in support for j in support if i!=j}),'exact_real_root_count_method':'SymPy Poly.count_roots(-oo, oo), Sturm sequence','exact_real_root_count':int(S.Poly(f,x).count_roots(-S.oo,S.oo)),'rational_factorization':[[str(q),int(e)] for q,e in S.factor_list(f)[1]],'discriminant':str(S.discriminant(f,x))}
 assert out['integer_binary_C'] and out['exact_real_root_count']==1
 Path('H7-a3b4-exact-feasibility.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
