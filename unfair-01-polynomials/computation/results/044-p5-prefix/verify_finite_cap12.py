#!/usr/bin/env python3
"""Fresh direct A*B convolution and isqrt floor; no degree-cap tail pruning."""
from pathlib import Path
from itertools import product
from math import isqrt
import hashlib,json,time
ROOT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials');HERE=Path(__file__).resolve().parent
INPUT=ROOT/'agents/coprime_saved_audit/round044-q2-search/revision2-corrected/x2-D12.json'
def plus(a,b):return(a[0]+b[0],a[1]+b[1])
def times(a,b):
 x=a[0]*b[0]+5*a[1]*b[1];y=a[0]*b[1]+a[1]*b[0];assert x%2==y%2==0;return(x//2,y//2)
def floor(z):
 u,v=z
 if v==0:assert u%2==0;return u//2
 q=isqrt(5*v*v)
 return (u+q)//2 if v>0 else (u-q-1)//2

def options(g,p,k):
 rest=(0,0)
 for j in range(1,min(k,5)+1):rest=plus(rest,(p[j][0]*g[k-j],p[j][1]*g[k-j]))
 f=floor(rest);choices=[-f] if rest[1] else [-f,1-f]
 return [(z,plus(rest,(2*z,0))) for z in choices]

def main():
 start=time.monotonic();doc=json.loads(INPUT.read_text());expected=list(product(range(-2,3),(-2,-1,1,2),range(-2,3),range(-2,3),range(-2,3)))
 observed=[tuple(r[k] for k in ['s','t','b','c','d']) for r in doc['rows']];assert observed==expected and len(expected)==2500
 all_rows=[];aggregate=[0]*12;last=0;tails=[]
 for row,key in zip(doc['rows'],expected):
  s,t,b,c,d=key;p=[(2,0),(2*b,0),(2*s+t,t),(2*c,0),(2*d,0),(2,0)];pb=[(u,-v) for u,v in p]
  states=[((1,),(1,),((2,0),),((2,0),))];layers=[]
  for k in range(1,13):
   fresh=[];tried=0
   for g,q,aa,bb in states:
    for ng,ak in options(g,p,k):
     for nq,bk in options(q,pb,k):
      tried+=1;az,bz=aa+(ak,),bb+(bk,);ck=(0,0)
      for i in range(k+1):ck=plus(ck,times(az[i],bz[k-i]))
      if ck in [(0,0),(2,0)]:fresh.append((g+(ng,),q+(nq,),az,bz))
   assert row['stages'][k-1]==[k,tried,len(fresh)],(key,k)
   layers.append(dict(k=k,tried=tried,accepted=[dict(G=list(g),Q=list(q)) for g,q,aa,bb in fresh]));states=fresh;aggregate[k-1]+=len(states)
   if states:last=max(last,k)
  for g,q,aa,bb in states:
   def factor(v,p):
    out=[(0,0)]*(len(v)+5)
    for i,x in enumerate(v):
     for j,z in enumerate(p):out[i+j]=plus(out[i+j],(x*z[0],x*z[1]))
    return out
   af,bf=factor(g,p),factor(q,pb)
   def unit(z):return floor(z)==0 if z[1] else z[0] in (0,2)
   failures=[dict(factor=nm,k=i,pair_over_two=list(z)) for nm,vv in [('A',af),('B',bf)] for i,z in enumerate(vv) if not unit(z)]
   cf=[(0,0)]*(len(af)+len(bf)-1)
   for i,x in enumerate(af):
    for j,y in enumerate(bf):cf[i+j]=plus(cf[i+j],times(x,y))
   cbad=[dict(k=i,pair_over_two=list(z)) for i,z in enumerate(cf) if z not in [(0,0),(2,0)]]
   assert failures or cbad,key
   tails.append(dict(parameters=list(key),G=list(g),Q=list(q),first_factor_failure=failures[0] if failures else None,first_product_failure=cbad[0] if cbad else None))
  all_rows.append(dict(parameters=list(key),free_prefix_layers=layers))
 out=dict(status='passed_complete_degree_cap12_search',mother_polynomials=2500,field='Q(sqrt5), omega=(1+sqrt5)/2',family='P=1+b*x+(s+t*omega)*x^2+c*x^3+d*x^4+x^5; s,b,c,d=-2..2;t=-2,-1,1,2',free_steps=12,forced_zero_tail_used=True,degree_cap_used=12,final_prefix_states=len(tails),all_final_prefixes_rejected=True,accepted_by_depth=aggregate,last_nonempty_depth=last,scope='Only these 2500 fixed P and integer constant-one G,Q with degrees at most 12; no arbitrary-degree or general exclusion.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),comparison_ledger_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),seconds=time.monotonic()-start)
 (HERE/'finite-cap12-summary.json').write_text(json.dumps(out,indent=2)+'\n');(HERE/'finite-cap12-all-states.json').write_text(json.dumps(all_rows,separators=(',',':'))+'\n');(HERE/'finite-cap12-tail-rejections.json').write_text(json.dumps(tails,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
