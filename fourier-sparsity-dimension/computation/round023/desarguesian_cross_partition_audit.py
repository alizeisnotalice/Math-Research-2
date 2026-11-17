#!/usr/bin/env python3
"""Exact finite-field audits for the R23 Desarguesian cross partition."""
import hashlib,json,platform,random,sys,time
from pathlib import Path
MOD={3:0b1011,4:0b10011,5:0b100101}
def mul(a,b,k):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a&(1<<k):a^=MOD[k]
 return out&((1<<k)-1)
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):x,y=a[j],a[j+h];a[j],a[j+h]=x+y,x-y
  h*=2
 return a
def rank(vs):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
def span(bs):
 s={0}
 for b in bs:s|={x^b for x in tuple(s)}
 return s
def audit(k):
 q=1<<k;n=2*k;N=1<<n;full=(1<<N)-1
 lines=[]
 for m in range(q):lines.append(sum(1<<(x|(mul(m,x,k)<<k)) for x in range(q)))
 lines.append(sum(1<<(y<<k) for y in range(q))) # vertical
 assert len(set(lines))==q+1 and all(x.bit_count()==q for x in lines)
 assert sum(lines)-(q+1)==full-1 # each line contains only the common origin
 neg=lines[:q//2];pos=lines[q//2:q];last=lines[q]
 f=0
 for L in neg:f^=L^1 # negative only off origin
 # General Walsh formula: E=sum eps=1, origin sign=1.
 W=fwht([1-2*((f>>x)&1) for x in range(N)]);assert set(abs(x) for x in W)=={q}
 blocks=[p^r for p,r in zip(neg,pos)]+[last]
 assert len(blocks)==q//2+1 and all(b.bit_count() in (2*q-2,q) for b in blocks)
 union=0
 for b in blocks:assert not(union&b);union|=b
 assert union==full and sum(blocks)==full
 coords=[sum(((x>>j)&1)<<x for x in range(N)) for j in range(n)]
 gens=[full]+coords+blocks;d=rank(gens);assert d==q//2+2*k+1
 assert rank(gens[:d])==d
 B=span(blocks);assert len(B)==1<<(q//2+1)
 # Every block choice preserves origin-sign=sum(line-signs); sample/full FWHT below.
 def check(codeword):return set(abs(x) for x in fwht([1-2*(((f^codeword)>>z)&1) for z in range(N)]))=={q}
 if k==3:
  V=span(gens);assert len(V)==1<<d and all(check(v) for v in V);mode='all full V phases'
 elif k==4:
  assert all(check(v) for v in B);assert all(check(v^g) for v in B for g in [0]+coords+[full]);mode='all block code plus every affine-basis translate'
 else:
  rng=random.Random(20261117)
  # samples are formed from actual generators, not arbitrary truth tables.
  samplewords=[]
  for _ in range(4097):
   c=rng.randrange(1<<d);w=0
   for j,g in enumerate(gens):
    if c>>j&1:w^=g
   samplewords.append(w)
  assert all(check(w) for w in samplewords);mode='4097 deterministic full-code samples'
 # A representative paired cross has two punctured spread spectra; record rank data.
 A=blocks[0];P,Q=neg[0],pos[0]
 PA=fwht([(1-2*((f>>x)&1)) if A>>x&1 else 0 for x in range(N)])
 supp={x for x,v in enumerate(PA) if v};assert len(supp)==2*(q-1) and set(PA)<={-q,0,q}
 return {'k':k,'q':q,'n':n,'field_modulus_bin':bin(MOD[k]),'spread_lines':q+1,'negative_lines':q//2,'positive_lines':q//2+1,'base_phase_bent':True,'block_count':len(blocks),'block_weights':[b.bit_count() for b in blocks],'blocks_partition_domain':True,'generator_rank':d,'expected_rank':q//2+2*k+1,'block_code_dimension':q//2+1,'verification_mode':mode,'seed':20261117 if k==5 else None,'phase_checks':(1<<d) if k==3 else (len(B)*(n+3) if k==4 else 4097),'representative_cross_hex':hex(A),'representative_cross_weight':A.bit_count(),'representative_partial_spectrum_size':len(supp),'representative_partial_values':sorted(set(PA))}
def main():
 t=time.monotonic();rows=[audit(k) for k in (3,4,5)]
 out={'schema':'a28-r023-desarguesian-cross-partition-v1','scope':'explicit Desarguesian spread construction k=3,4,5; formula proves all block flips, finite FWHT scope reported separately','arithmetic':'exact GF(2) field multiplication and integer FWHT','rows':rows,'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.monotonic()-t,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 Path(__file__).with_name('desarguesian-cross-partition-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
