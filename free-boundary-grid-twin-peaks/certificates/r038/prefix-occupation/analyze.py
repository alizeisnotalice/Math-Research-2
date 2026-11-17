#!/usr/bin/env python3
"""Exact diagnostics of a saved finite one-peak word, never a probability sample."""
import argparse,json,hashlib,pathlib,math
from fractions import Fraction
p=argparse.ArgumentParser();p.add_argument('word',type=pathlib.Path);p.add_argument('--output',type=pathlib.Path,required=True);z=p.parse_args();lines=z.word.read_text().splitlines();h=dict(zip(lines[0].split()[::2],lines[0].split()[1::2]));n=int(h['n']);N=n*n;r=int(h['root']);w=list(map(int,lines[1].split()[1:]));assert len(w)==N and set(w)==set(range(N)) and w[0]==r
pos=[0]*N
for i,v in enumerate(w):pos[v]=i+1
def nb(v):
 x,y=v%n,v//n
 return [Y*n+X for X,Y in [(x-1,y),(x+1,y),(x,y-1),(x,y+1)] if 0<=X<n and 0<=Y<n]
early=[[u for u in nb(v) if pos[u]<pos[v]] for v in range(N)];assert all(early[v] for v in range(N) if v!=r)
terms=[];psi=Fraction();geo=Fraction();far={4:Fraction(),2:Fraction(),1:Fraction()};sentinel=Fraction();layers={};manhattan=Fraction();far_e={k:Fraction() for k in far};far_m={k:Fraction() for k in far}
for v in range(N):
 if v==r or v in nb(r):continue
 l=min(pos[u] for u in nb(v));b=min((pos[u] for u in nb(v) if early[u]==[v]),default=N+1);R=b-l-1;assert R>0
 t=Fraction(l-1,R);d=abs(v%n-r%n)+abs(v//n-r//n);e=Fraction((v%n-r%n)**2+(v//n-r//n)**2,N)
 psi+=t;geo+=t*e;manhattan+=t*Fraction(d*d,N)
 if b==N+1:sentinel+=t
 for den in far:
  if den*d>n:
   far[den]+=t;far_e[den]+=t*e;far_m[den]+=t*Fraction(d*d,N)
 block=(l-1).bit_length();layers[block]=layers.get(block,Fraction())+t
 terms.append(dict(vertex=v,rank=pos[v],l=l,b=b,R=R,distance=d,term=str(t)))
assert abs(float(psi)-float(h['psi']))<=1e-11*max(1,float(psi))
out=dict(scope='finite deterministic candidate only',input_sha256=hashlib.sha256(z.word.read_bytes()).hexdigest(),header=h,n=n,root=r,psi=str(psi),psi_float=float(psi),normalized=float(psi)/(N*math.log(n)),psi_euclidean=str(geo),psi_manhattan_squared=str(manhattan),far_psi_euclidean={str(k):str(v) for k,v in far_e.items()},far_psi_manhattan_squared={str(k):str(v) for k,v in far_m.items()},far_psi={str(k):str(v) for k,v in far.items()},sentinel_psi=str(sentinel),rank_dyadic_layers={str(k):str(v) for k,v in sorted(layers.items())},terms=terms)
hist={}
for x in terms:hist[x['l']]=hist.get(x['l'],Fraction())+Fraction(1,x['R'])
g=sum(hist.values(),Fraction());occupation=[]
for t in range(1,N):
 occupation.append(g);g-=hist.get(t+1,Fraction())
assert sum(occupation,Fraction())==psi and all(occupation[i]>=occupation[i+1] for i in range(len(occupation)-1))
out['prefix_occupation']={'gamma':str(occupation[0]),'sum':str(sum(occupation,Fraction())),'first_dominating_prefix':next((i+1 for i,v in enumerate(occupation) if v==0),None),'identity_pass':True}
z.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['terms','header']},indent=2))
