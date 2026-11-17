#!/usr/bin/env python3
"""Complete exact audit of 2q switches from f=x0x1+x2x3 on F2^4."""
import hashlib,json,platform,time
from pathlib import Path

OUT=Path(__file__).with_name('n4-double-switch-audit.json')
N=16
def dot(a,b):return (a&b).bit_count()&1
def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
def bent(mask):return all(abs(x)==4 for x in fwht([1-2*((mask>>i)&1) for i in range(N)]))
def rank(xs):
 p={}
 for x in xs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
def span(bs):
 o={0}
 for b in bs:o|={x^b for x in tuple(o)}
 return o
def affine_flats(dim):
 seen=set();out=[]
 for basis_code in range(1<<(4*dim)):
  bs=[(basis_code>>(4*j))&15 for j in range(dim)]
  if rank(bs)!=dim:continue
  H=frozenset(span(bs))
  for p in range(16):
   S=frozenset(p^x for x in H)
   if S not in seen:seen.add(S);out.append(S)
 return out
def maskof(S):return sum(1<<x for x in S)
def is_affine(S):
 p=next(iter(S)); H={x^p for x in S}; return len(S)==(1<<rank(H)) and span(list(H))==H
def main():
 start=time.monotonic();f0=sum((((x&1)*((x>>1)&1))^(((x>>2)&1)*((x>>3)&1)))<<x for x in range(16))
 assert bent(f0)
 flats2=affine_flats(2);assert len(flats2)==140
 legal2=[maskof(S) for S in flats2 if bent(f0^maskof(S))]
 assert len(legal2)>0
 bends=[g for g in range(1<<16) if bent(g)]
 w8=[];records=[]
 for g in bends:
  v=f0^g
  if v.bit_count()!=8:continue
  contained=[b for b in legal2 if b&~v==0]
  decomps=[]
  for b in contained:
   c=v^b
   if c in legal2 and b<c:decomps.append((b,c))
  rec={'switch_mask':v,'support_affine':is_affine({x for x in range(16) if v>>x&1}),
       'contained_legal_q_flats':len(contained),'decomposition_count':len(decomps)}
  if decomps:rec['first_decomposition']=list(decomps[0])
  records.append(rec);w8.append(v)
 weight_hist={}
 for g in bends:weight_hist[str((f0^g).bit_count())]=weight_hist.get(str((f0^g).bit_count()),0)+1
 non=[r for r in records if not r['decomposition_count']]
 result={'schema':'a28-r018-n4-double-switch-v1','base_mask':f0,'total_bent_phases':len(bends),
  'weight_histogram':weight_hist,'affine_q_flats_total':len(flats2),'legal_q_flats_total':len(legal2),
  'weight_2q_switches':len(records),'weight_2q_affine_supports':sum(r['support_affine'] for r in records),
  'weight_2q_decomposable':len(records)-len(non),'weight_2q_nondecomposable':len(non),
  'nondecomposable_certificates':non[:10],'elapsed_seconds':time.monotonic()-start,
  'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':platform.python_version(),'arithmetic':'integer FWHT and exhaustive affine-flat scan','failed_assertions':0}
 OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
