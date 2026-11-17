#!/usr/bin/env python3
"""Independent complete box/partition replay; descending integer disk recursion."""
import hashlib,json,math,time,platform
from pathlib import Path
from collections import Counter
ROOT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/coprime_stretches/round045-core5')
HERE=Path(__file__).resolve().parent
POINTS=((1,4),(1,2),(3,4),(1,1),(4,3),(2,1),(4,1))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def hpoly(t,n,c,d):
 f=[1,0,0,c,d,1];h=[0]*11
 for i,u in enumerate(f):
  for j,v in enumerate(f):h[i+j]+=u*v
 for i,u in enumerate(f):h[i+2]+=t*u
 h[4]+=n
 return h
def disk_failure(ascending):
 # Work in descending order, cancel the last coefficient and remove it.
 a=list(reversed(ascending))
 while len(a)>1:
  content=math.gcd(*a)
  if a[0]<0:content=-content
  a=[x//content for x in a]
  if abs(a[-1])>=a[0]:return (len(a)-1,abs(a[-1]),a[0])
  lead,last=a[0],a[-1]
  a=[lead*a[i]-last*a[-1-i] for i in range(len(a)-1)]
 return None
def axis_failure(t,n,c,d):
 delta=t*t-4*n
 for p,q in POINTS:
  # Exact Horner evaluation of 2P's rational part after clearing q^5.
  coeff=[2,0,t,2*c,2*d,2]
  L=sum(v*p**i*q**(5-i) for i,v in enumerate(coeff))
  margin=L*L-delta*p**4*q**6
  if L<=0 or margin<=0:return {'point':[p,q],'L':L,'product_numerator_x4':margin}
 return None
def main():
 started=time.monotonic();summary=json.loads((ROOT/'core5_summary.json').read_text())
 for name,key in [('core5_survivors.json','survivors_sha256'),('core5_exclusion_witnesses.jsonl','witnesses_sha256'),('core5_delta_family_rejections.jsonl','delta_family_witnesses_sha256')]:
  assert sha(ROOT/name)==summary[key]
 kept=json.loads((ROOT/'core5_survivors.json').read_text());keep_set={tuple(x) for x in kept};assert len(keep_set)==len(kept)
 families=iter((ROOT/'core5_delta_family_rejections.jsonl').read_text().splitlines())
 counts=Counter();survivors=[]
 with (ROOT/'core5_exclusion_witnesses.jsonl').open() as records:
  for t in range(-19,20):
   for n in range(-99,100):
    delta=t*t-4*n;counts['raw_box']+=513
    bad=delta<=0 or math.isqrt(delta)**2==delta or delta>=(20-abs(t))**2
    if bad:
     e=json.loads(next(families));assert e=={'T':t,'N':n,'stage':'delta_trace','delta':delta,'threshold':(20-abs(t))**2,'covered_c_d_tuples':513}
     continue
    for c in range(-13,14):
     for d in range(-9,10):
      row=(t,n,c,d);counts['delta_trace']+=1;w=axis_failure(t,n,c,d)
      if w is not None:stage='positive_axis'
      else:
       counts['positive_axis']+=1;h=hpoly(*row);stage='annulus_schur'
       for label,coef in [('H(2x)',h),('reverse(H)(2x)',h[::-1])]:
        ans=disk_failure([v*2**i for i,v in enumerate(coef)])
        if ans is not None:
         w=dict(orientation=label,degree=ans[0],abs_constant=ans[1],abs_leading=ans[2]);break
      if w is None:
       assert row in keep_set;counts['annulus_schur']+=1;survivors.append(list(row))
      else:
       e=json.loads(next(records));assert e=={'candidate':list(row),'stage':stage,'witness':w},row
  assert not records.read().strip()
 assert next(families,None) is None and survivors==kept and dict(counts)==summary['counts']
 out=dict(status='passed_full_independent_integer_partition',counts=dict(counts),python=platform.python_version(),source_sha256=sha(Path(__file__)),survivors_sha256=sha(ROOT/'core5_survivors.json'),summary_sha256=sha(ROOT/'core5_summary.json'),seconds=time.monotonic()-started,scope='404 necessary P5 norm-core tuples; arbitrary-degree cofactors remain unresolved.')
 (HERE/'parent-summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
