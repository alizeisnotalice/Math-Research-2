#!/usr/bin/env python3
"""Independent exact checker for any recorded affine-flat witnesses."""
import argparse,json,hashlib
from pathlib import Path
def rank(vs,n):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
def mask(n,basis,off):
 ans=0
 for x in range(1<<n):
  if all(((x&a).bit_count()&1)==((off>>j)&1) for j,a in enumerate(basis)):ans|=1<<x
 return ans
ap=argparse.ArgumentParser();ap.add_argument('result',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();raw=a.result.read_bytes();d=json.loads(raw)
checks=[]
for row in d['bounded_k3_searches']:
 flats=row.get('flats',row.get('best_flats',[])); n,k=row['n'],row['k']; used=0; norms=[]
 for f in flats:
  b=f['normal_basis'];m=mask(n,b,f['offset_bits'])
  if m!=int(f['point_mask_hex'],16) or m.bit_count()!=1<<(n-k) or used&m:raise AssertionError('flat')
  used|=m;norms+=b
 if len(flats)>1<<k:raise AssertionError('m')
 checks.append({'n':n,'k':k,'m':len(flats),'normal_rank':rank(norms,n)})
a.out.write_text(json.dumps({'status':'passed','input_sha256':hashlib.sha256(raw).hexdigest(),'witness_checks':checks},indent=2)+'\n')
