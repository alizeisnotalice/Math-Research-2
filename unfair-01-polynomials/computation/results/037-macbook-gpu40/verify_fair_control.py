"""Recompute the fair control with exact integer convolution."""
from pathlib import Path
import json,hashlib
path=Path(__file__).with_name('next48-fair-control.json');data=json.loads(path.read_text())
A=[round(x) for x in data['sample']['A']];B=[round(x) for x in data['sample']['B']]
assert set(A)|set(B)<={0,1} and A[0]==B[0]==1
assert [i for i,v in enumerate(A) if v]==[0,1,3]
assert [i for i,v in enumerate(B) if v]==[0,5,11,20,45]
C=[0]*(len(A)+len(B)-1)
for i,a in enumerate(A):
    for j,b in enumerate(B):C[i+j]+=a*b
assert set(C)<={0,1}
assert [i for i,v in enumerate(C) if v]==[0,1,3,5,6,8,11,12,14,20,21,23,45,46,48]
assert [i for i,v in enumerate(C) if v]==data['C_support']
err=max(abs(x-y) for v,w in [(data['sample']['A'],A),(data['sample']['B'],B)] for x,y in zip(v,w))
assert err<1e-5 and err==data['maximum_rounding_change']
print(json.dumps({'status':'passed','binary_factors':True,'binary_product':True,'maximum_rounding_change':err,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'scope':'exact fair control only; not an unfair counterexample'}))
