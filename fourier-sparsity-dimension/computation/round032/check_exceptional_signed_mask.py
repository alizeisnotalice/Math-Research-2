"""Exhaust signs on the exceptional support modulo affine characters.

Eight affine-independent points are fixed positive. W(0) then restricts
the remaining 24 negative signs to weights 8, 16 or 24. All 128 Walsh
coefficients are checked by exact popcount, with short-circuit rejection.
"""
import itertools,json,time,hashlib,sys,platform
from pathlib import Path
P=Path(__file__).resolve().parent
S=[x for x in range(128) if (((x>>0&1)*(x>>1&1)*(x>>3&1))^((x>>0&1)*(x>>4&1)*(x>>5&1))^((x>>1&1)*(x>>2&1)*(x>>6&1)))]
assert len(S)==32
basis={};fixed=[]
for i,x in enumerate(S):
 v=x|128
 while v:
  p=v.bit_length()-1
  if p in basis:v^=basis[p]
  else:basis[p]=v;fixed.append(i);break
assert len(fixed)==8
free=[i for i in range(32) if i not in fixed];F=sum(1<<i for i in free)
C=[sum(1<<i for i,x in enumerate(S) if (a&x).bit_count()%2) for a in range(128)]
C.sort(key=lambda x:abs(x.bit_count()-16),reverse=True)
t=time.time();checked=0;solutions=[]
for inds in itertools.combinations(free,8):
 n=sum(1<<i for i in inds)
 for N in (n,F^n):
  checked+=1
  if all((N^c).bit_count() in (8,16,24) for c in C):solutions.append(N)
checked+=1
if all((F^c).bit_count() in (8,16,24) for c in C):solutions.append(F)
assert checked==1470943
out=dict(python=sys.version,platform=platform.platform(),points=S,fixed_positive_indices=fixed,free_indices=free,checked=checked,solutions=solutions,seconds=time.time()-t,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(P/'exceptional_signed_mask_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
