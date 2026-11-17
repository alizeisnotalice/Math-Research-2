#!/usr/bin/env python3
"""Independent explicit two-branch prefix audit for degree6 irrational position1."""
import math,json,gzip,time,hashlib
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent
def add(x,y):return x[0]+y[0],x[1]+y[1]
def scale(a,x):return a*x[0],a*x[1]
def mul(x,y,t,n):return x[0]*y[0]-n*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+t*x[1]*y[1]
def options(x,t,delta):
 p,q=x
 if not q:return [(-p,(0,0)),(1-p,(1,0))]
 u=2*p+q*t;s=math.isqrt(q*q*delta)
 f=(u+s)//2 if q>0 else (u-s-1)//2
 return [(-f,(p-f,q))]
def convolution(a,b,k,t,n):
 s=(0,0)
 for i in range(k+1):s=add(s,mul(a[i],b[k-i],t,n))
 return s
def bit(x):return x in ((0,0),(1,0))
def main():
 start=time.monotonic();counts=Counter();survivors=[];pairs=[]
 with gzip.open(HERE/'single-terminal-witnesses.jsonl.gz','wt') as journal:
  for t in range(-14,15):
   for n in range(-56,57):
    delta=t*t-4*n
    if delta<=0 or math.isqrt(delta)**2==delta or delta>=(15-abs(t))**2:continue
    pairs.append([t,n]);counts['tn_pairs']+=1;a=(0,1);b=(t,-1)
    (g1,A1),=options(a,t,delta);(q1,B1),=options(b,t,delta)
    assert bit(add(A1,B1));counts['pass_c1']+=1
    for c in range(-21,22):
     stage2=[]
     for g2,A2 in options(add((c,0),scale(g1,a)),t,delta):
      for q2,B2 in options(add((c,0),scale(q1,b)),t,delta):
       if bit(add(add(A2,B2),mul(A1,B1,t,n))):stage2.append((g2,q2,A2,B2))
     if not stage2:
      journal.write(json.dumps([t,n,c,'C2_all_choices_fail'],separators=(',',':'))+'\n');continue
     assert len(stage2)==1;counts['pass_c2']+=1
     g2,q2,A2,B2=stage2[0]
     for d in range(-30,31):
      stage3=[]
      for g3,A3 in options(add((d+c*g1,0),scale(g2,a)),t,delta):
       for q3,B3 in options(add((d+c*q1,0),scale(q2,b)),t,delta):
        AA=[(1,0),A1,A2,A3];BB=[(1,0),B1,B2,B3]
        if bit(convolution(AA,BB,3,t,n)):stage3.append((g3,q3,AA,BB))
      if not stage3:
       journal.write(json.dumps([t,n,c,d,'C3_all_choices_fail'],separators=(',',':'))+'\n');continue
      assert len(stage3)==1;counts['pass_c3']+=1
      g3,q3,AA,BB=stage3[0]
      ka=add((c*g2+d*g1,0),scale(g3,a));kb=add((c*q2+d*q1,0),scale(q3,b))
      ea=dict(options(ka,t,delta));eb=dict(options(kb,t,delta))
      common=sorted(set(ea)&set(eb)&set(range(-21,22)));checks=[]
      for e in common:
       C4=convolution(AA+[ea[e]],BB+[eb[e]],4,t,n);checks.append([e,C4])
       if bit(C4):survivors.append([t,n,c,d,e,g1,g2,g3,q1,q2,q3])
      journal.write(json.dumps([t,n,c,d,'C4_endpoint_G4_Q4_zero',ka,kb,sorted(ea),sorted(eb),checks],separators=(',',':'))+'\n')
 assert len(pairs)==455 and not survivors
 source=HERE.parent/'agents/algebraic_global/round046/core6_single/core6_single_results.json'
 other=json.loads(source.read_text())
 for k in ('pass_c1','pass_c2','pass_c3'):assert counts[k]*2==other['counts'][k]
 out={'status':'passed_complete_single_core_prefix_exclusion','counts':dict(counts),'survivors':survivors,'tn_pairs':pairs,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'witnesses_sha256':hashlib.sha256((HERE/'single-terminal-witnesses.jsonl.gz').read_bytes()).hexdigest(),'seconds':time.monotonic()-start,'scope':'All degree6 single-irrational-position1 cores under the proved root envelope; arbitrary-degree cofactors covered by necessary G4=Q4=0. No linked-position2,3 claim.'}
 (HERE/'single-core-summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='tn_pairs'}))
if __name__=='__main__':main()
