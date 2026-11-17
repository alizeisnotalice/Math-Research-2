#!/usr/bin/env python3
"""Direct AB exact prefix tree for linked six-degree cores; endpoint rows 1,3=0."""
import json,math,hashlib,time,platform,gzip
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent;SOURCE=HERE/'linked-prefix-input.json';LIMIT=128
def mul(x,y,t,n):
 p,q=x;r,s=y
 return p*r-n*q*s,p*s+q*r+t*q*s
def choices(p,q,t,delta):
 if not q:return [(-p,(0,0)),(1-p,(1,0))]
 u=2*p+q*t;s=math.isqrt(q*q*delta);f=(u+s)//2 if q>0 else (u-s-1)//2
 return [(-f,(p-f,q))]
def main():
 start=time.monotonic();cores=json.loads(SOURCE.read_text())['rows'];rows=[];counts=Counter()
 with gzip.open(HERE/'linked-direct-trees.jsonl.gz','wt') as out:
  for index,core in enumerate(cores):
   b,t,n,ell,d,e=core;delta=t*t-4*n
   P=((1,0),(b,0),(0,1),(ell,b),(d,0),(e,0),(1,0))
   Qp=tuple((p+t*q,-q) for p,q in P)
   states=[((1,),(1,),((1,0),),((1,0),))];levels=[];tree=[];status='unknown_depth'
   for k in range(1,LIMIT+1):
    nxt=[];parents=[]
    for G,Q,A,B in states:
     ap=aq=bp=bq=0
     for j in range(1,min(k,6)+1):
      ap+=P[j][0]*G[k-j];aq+=P[j][1]*G[k-j];bp+=Qp[j][0]*Q[k-j];bq+=Qp[j][1]*Q[k-j]
     inner=[0,0]
     for j in range(1,k):
      u,v=mul(A[j],B[k-j],t,n);inner[0]+=u;inner[1]+=v
     accepted=[];rejected=[]
     for g,a in choices(ap,aq,t,delta):
      for q,v in choices(bp,bq,t,delta):
       c=(inner[0]+a[0]+v[0],inner[1]+a[1]+v[1]);assert c[1]==0
       if k in (1,3) and (a!=(0,0) or v!=(0,0)):
        rejected.append([g,q,'endpoint_zero',a,v]);continue
       if c[0] not in (0,1):rejected.append([g,q,'C_nonbinary',c[0]]);continue
       accepted.append([g,q]);nxt.append((G+(g,),Q+(q,),A+(a,),B+(v,)))
     parents.append({'G':G,'Q':Q,'accepted':accepted,'rejected':rejected})
    assert len(nxt)==len({(x[0],x[1]) for x in nxt})
    levels.append(len(nxt));tree.append({'k':k,'parents':parents});states=nxt
    if not states:status='dead_all_prefixes';break
    if len(states)>100000:status='unknown_state_cap';break
   row={'core':core,'status':status,'depth':len(levels),'counts':levels}
   if states:row['surviving_G_Q']=[[G,Q] for G,Q,A,B in states]
   rows.append(row);counts[status]+=1;out.write(json.dumps({'core':core,'levels':tree},separators=(',',':'))+'\n')
   if (index+1)%1000==0:print({'done':index+1,'statuses':dict(counts)},flush=True)
 summary={'status_counts':dict(counts),'total':len(cores),'rows':rows,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,'python':platform.python_version(),'scope':'Requires proved root-input coverage and endpoint A1=B1=A3=B3=0; empty prefix trees exclude any finite cofactor degrees. Unknowns remain open.'}
 (HERE/'linked-direct-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)
if __name__=='__main__':main()
