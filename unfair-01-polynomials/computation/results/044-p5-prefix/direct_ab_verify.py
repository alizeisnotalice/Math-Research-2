#!/usr/bin/env python3
"""Independent direct-A*B semantic verifier for corrected x2 prefix ledgers."""
import argparse,json,hashlib,time
from pathlib import Path
def s(u,v):
 if not u and not v:return 0
 if u>=0 and v>=0:return 1
 if u<=0 and v<=0:return -1
 return 1 if (u*u>5*v*v if u>0 else 5*v*v>u*u) else -1
def c(z,n):return s(z[0]-2*n,z[1])
def add(*z):return(sum(x for x,y in z),sum(y for x,y in z))
def I(n):return(2*n,0)
def M(z,n):return z[0]*n,z[1]*n
def box(z):return c(z,0)>=0 and c(z,1)<=0
def floor(z):
 r=(abs(z[0])+6*abs(z[1]))//2+2;l,h=-r,r
 while l<h:
  m=(l+h+1)//2
  if c(z,m)>=0:l=m
  else:h=m-1
 return l
def coeff(g,p,k):return add(*(M(p[j],g[k-j]) for j in range(min(k,5)+1) if k-j<len(g)))
def opt(g,p):
 k=len(g); rest=add(*(M(p[j],g[k-j]) for j in range(1,min(k,5)+1)))
 f=floor((-rest[0],-rest[1]));return [x for x in(f,f+1) if box(add(I(x),rest))]
def directC(g,q,p,pb,k):
 return sum(1 for i in range(k+1) for _ in [0] if box(coeff(g,p,i)) and box(coeff(q,pb,k-i)) and coeff(g,p,i)!=(0,0) and coeff(q,pb,k-i)!=(0,0)) # replaced below
def run(row,D):
 s0,t,b,cc,d=(row[x] for x in('s','t','b','c','d'));a=(2*s0+t,t);ab=(2*s0+t,-t);p=(I(1),I(b),a,I(cc),I(d),I(1));pb=(I(1),I(b),ab,I(cc),I(d),I(1));S=[((1,),(1,))];L=[]
 for k in range(1,D+1):
  N=[];pre=0
  for g,q in S:
   for x in opt(g,p):
    for y in opt(q,pb):
     pre+=1;gg,qq=g+(x,),q+(y,); A=[coeff(gg,p,i) for i in range(k+1)];B=[coeff(qq,pb,i) for i in range(k+1)]
     # Direct K multiplication: product pair is rational here; assert exactly 0 or 1.
     z=(0,0)
     for i in range(k+1):
      u,v=A[i];r,w=B[k-i];z=add(z,((u*r+5*v*w)//2,(u*w+v*r)//2))
     if z in ((0,0),(2,0)):N.append((gg,qq))
  L.append([k,pre,len(N)]);S=N
 return L,len(S)
ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--output',type=Path,required=True);ns=ap.parse_args();t=time.monotonic();doc=json.loads(ns.input.read_text());err=[]
if doc.get('processed')!=2500 or len(doc.get('rows',[]))!=2500:err.append({'kind':'coverage'})
for i,r in enumerate(doc.get('rows',[])):
 L,n=run(r,doc['cap'])
 if r.get('stages',[])[:doc['cap']]!=L:err.append({'kind':'stage','row':i})
 if r.get('free_prefix_states_after_cap')!=n:err.append({'kind':'free','row':i})
out={'status':'passed' if not err else 'rejected','errors':err,'rows':len(doc.get('rows',[])),'input_sha256':hashlib.sha256(ns.input.read_bytes()).hexdigest(),'seconds':time.monotonic()-t}
ns.output.write_text(json.dumps(out,indent=2)+'\n')
if err:raise SystemExit(2)
