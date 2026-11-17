#!/usr/bin/env python3
"""Enumerate all 115 unital algebras in the fixed D7 code and build completed selectors.
The returned independent basis is in echelon form; global code deduplication is not claimed.
"""
import argparse,itertools,json
from pathlib import Path
parser=argparse.ArgumentParser();parser.add_argument("--out",type=Path,required=True);args=parser.parse_args()
def sp(a):
 o={0}
 for x in a:o|={z^x for z in list(o)}
 return o
def rk(a):
 b={}
 for x in a:
  while x:
   q=x.bit_length()-1
   if q in b:x^=b[q]
   else:b[q]=x;break
 return len(b)
B=[]
for i in range(7):
 m=0
 for p in range(16):
  z0,z1,w0,w1=[p>>j&1 for j in range(4)]
  if [1,w0,w1,w0*w1,z0,z1,z1*w0][i]:m|=1<<p
 B.append(m)
def ev(u):
 r=0
 for i in range(7):
  if u>>i&1:r^=B[i]
 return r
L={ev(u):u for u in range(128)}
def mul(a,b):return L.get(ev(a)&ev(b))
subs=[]
for d in range(8):
 for piv in itertools.combinations(range(7),d):
  free=[(i,j)for i,x in enumerate(piv)for j in range(x+1,7)if j not in piv]
  for q in range(1<<len(free)):
   rr=[1<<x for x in piv]
   for k,(i,j) in enumerate(free):
    if q>>k&1:rr[i]|=1<<j
   E=sp(rr)
   if 1 not in E or any(mul(a,b) not in E for a in E for b in E):continue
   core=[u for u in range(128)if all(mul(e,u)!=None for e in E)]
   G=[]
   for w in range(2):
    for u in core:G.append(sum(((ev(u)>>p)&1)<<(p|(z<<4)|(w<<5))for p in range(16)for z in range(2)))
   for e in E:G.append(sum(((ev(e)>>p)&1)<<(p|16|(w<<5))for p in range(16)for w in range(2)))
   for j in range(4):G.append(sum(((p>>j)&1)<<(p|z<<4|w<<5)for p in range(16)for z in range(2)for w in range(2)))
   I=[]
   for g in G:
    if rk(I+[g])>len(I):I.append(g)
   subs.append((rr,tuple(sorted(I))))
codes={}
for e,c in subs:codes.setdefault(c,[]).append(e)

def echelon_basis(v):
 b={}
 for x in v:
  while x:
   q=x.bit_length()-1
   if q in b:x^=b[q]
   else:b[q]=x;break
 return [b[k] for k in sorted(b,reverse=True)]
json.dump({"schema":"r014-selector-affine-v1","domain_bits":"x0,x1,y0,y1,z_new,w_new","records":[{"rank":len(c),"basis_masks":echelon_basis(c),"E_rows":es} for c,es in codes.items()]},args.out.open("w"),indent=2)

assert len(codes)==115
print("generated",len(codes))
