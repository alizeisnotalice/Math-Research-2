#!/usr/bin/env python3
"""Exact direct A*B prefix tree, independently of H*G*Q search.

All finite solutions extend by zero and hence occur in this necessary prefix
tree. Exhaustion certifies arbitrary-degree absence for that particular core;
surviving prefixes certify neither finite solvability nor absence.
"""
import json,math,hashlib,time,platform
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/coprime_stretches/round045-core5/core5_survivors.json')
LIMIT=64
def floor_pair(p,q,t,delta):
 if not q:return p
 u=2*p+q*t;s=math.isqrt(q*q*delta)
 return (u+s)//2 if q>0 else (u-s-1)//2
def choices(p,q,t,delta):
 if not q:return [(-p,(0,0)),(1-p,(1,0))]
 z=-floor_pair(p,q,t,delta);return [(z,(p+z,q))]
def mul(x,y,t,n):
 p,q=x;r,s=y
 return p*r-n*q*s,p*s+q*r+t*q*s
def extend(state,k,t,n,c,d):
 G,Q,A,B=state;delta=t*t-4*n
 get=lambda f,i:f[i] if i>=0 else 0
 ga=get(G,k-2);qb=get(Q,k-2)
 ag=c*get(G,k-3)+d*get(G,k-4)+get(G,k-5)
 bg=t*qb+c*get(Q,k-3)+d*get(Q,k-4)+get(Q,k-5)
 inner=[0,0]
 for i in range(1,k):
  u,v=mul(A[i],B[k-i],t,n);inner[0]+=u;inner[1]+=v
 out=[];bad=[]
 for g,a in choices(ag,ga,t,delta):
  for q,b in choices(bg,-qb,t,delta):
   coeff=(inner[0]+a[0]+b[0],inner[1]+a[1]+b[1])
   assert coeff[1]==0,(k,G,Q,coeff)
   if coeff[0] in (0,1):out.append((G+(g,),Q+(q,),A+(a,),B+(b,)))
   else:bad.append({'g':g,'q':q,'product_coefficient':coeff[0]})
 return out,bad
def main():
 start=time.monotonic();cores=json.loads(SOURCE.read_text());rows=[];trees=[]
 for index,(t,n,c,d) in enumerate(cores):
  states=[((1,),(1,),((1,0),),((1,0),))];counts=[];tree=[];status='alive_at_depth_limit'
  for k in range(1,LIMIT+1):
   nxt=[];stage=[]
   for state in states:
    good,bad=extend(state,k,t,n,c,d);nxt.extend(good)
    stage.append({'G':list(state[0]),'Q':list(state[1]),'accepted_last':[[s[0][-1],s[1][-1]] for s in good],'rejected':bad})
   assert len({(s[0],s[1]) for s in nxt})==len(nxt)
   states=nxt;counts.append(len(states));tree.append({'k':k,'parents':stage})
   if not states:status='all_prefixes_exhausted';break
   if len(states)>10000:status='state_cap_unknown';break
  row={'core':[t,n,c,d],'status':status,'depth':len(counts),'counts':counts,'surviving_G_Q':[[list(s[0]),list(s[1])] for s in states]}
  rows.append(row);trees.append({'core':row['core'],'tree':tree})
  if (index+1)%100==0:print(json.dumps({'cores_done':index+1,'dead':sum(r['status']=='all_prefixes_exhausted' for r in rows)}),flush=True)
 out={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':platform.python_version(),'depth_limit':LIMIT,'rows':rows,'dead':sum(r['status']=='all_prefixes_exhausted' for r in rows),'total':len(cores),'seconds':time.monotonic()-start}
 (HERE/'direct-prefix-summary.json').write_text(json.dumps(out,indent=2)+'\n')
 import gzip
 (HERE/'direct-prefix-trees.json.gz').write_bytes(gzip.compress(json.dumps(trees,separators=(',',':')).encode(),mtime=0))
 print(json.dumps({k:v for k,v in out.items() if k!='rows'}),flush=True)
if __name__=='__main__':main()
