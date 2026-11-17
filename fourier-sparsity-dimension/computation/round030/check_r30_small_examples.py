#!/usr/bin/env python3
"""Exact integer-FWHT checks for the two explicit R30 nonparallel examples."""
from __future__ import annotations
import hashlib,json,platform,sys,time
from pathlib import Path
OUT=Path(__file__).with_name("r30_small_examples_result.json")
Q,q=64,8
def par(x):return x.bit_count()&1
def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for s in range(0,len(a),2*h):
   for j in range(h):a[s+j],a[s+j+h]=a[s+j]+a[s+j+h],a[s+j]-a[s+j+h]
  h*=2
 return a
def rank(vs):
 b={}
 for v in vs:
  while v:
   p=v.bit_length()-1
   if p in b:v^=b[p]
   else:b[p]=v;break
 return len(b)
def f(x):
 a,b,c,d,u,v=((x>>i)&1 for i in range(6));return (a&b)^(c&d)^(u&v)
def oldmask(mask,x):
 a,b,c,d=((x>>i)&1 for i in range(4));vals=(1,a,b,c,d,a&c,a&d)
 return par(sum(vals[i]<<i for i in range(7))&mask)
def gen(mask,x):return oldmask(mask,x) if not((x>>5)&1) else 0
def setmask(pred):return sum(1<<x for x in range(Q) if pred(x))
def independent(vs):
 ans=[];last=0
 for v in vs:
  if rank(ans+[v])>last:ans.append(v);last+=1
 return ans
A0=setmask(lambda x:not(x&1) and not((x>>3)&1) and not((x>>5)&1))
A1=setmask(lambda x:not(x&1) and ((x>>3)&1) and not((x>>5)&1))
A2=setmask(lambda x:(x&1) and not((x>>2)&1) and not((x>>5)&1))
A3=setmask(lambda x:(x&1) and ((x>>2)&1) and not((x>>5)&1))
def flat(block):
 pts=[x for x in range(Q) if block>>x&1];base=pts[0];directions=independent([x^base for x in pts])
 span={0}
 for d in directions:span|={x^d for x in list(span)}
 return base,directions,span
def span(gens):
 s={0}
 for g in gens:s|={x^g for x in list(s)}
 return s
def verify(name,genmasks,blocks,blockcoords,expectedrank,expectedslope,expected_size,require_full_D):
 words=[setmask(lambda x,m=m:gen(m,x)) for m in genmasks]
 assert rank(words)==expectedrank
 flats=[flat(b) for b in blocks]
 assert all(len(pts)==8 and block.bit_count()==8 for (_,_,pts),block in zip(flats,blocks))
 assert all({base^x for x in pts}=={x for x in range(Q) if block>>x&1} for (base,_,pts),block in zip(flats,blocks))
 shapes=[frozenset(pts) for _,_,pts in flats]
 assert len(set(shapes))==2
 R=0
 for block in blocks:
  assert not R&block
  R|=block
 assert all(not (word & (((1<<Q)-1)^R)) for word in words)
 rebuilt=[]
 for coord in blockcoords:
  x=0
  for j,w in enumerate(words):
   if coord>>j&1:x^=w
  rebuilt.append(x)
 assert rebuilt==blocks
 labels=[];fixed=None;directions=None;allbent=partialok=True
 for p in range(1<<len(words)):
  signs=[]
  for x in range(Q):
   bit=f(x)
   for j,m in enumerate(genmasks):
    if p>>j&1:bit^=gen(m,x)
   signs.append(1 if not bit else -1)
  allbent &= all(abs(v)==q for v in fwht(signs))
  coverage=[0]*Q;dnow=[]
  for block in blocks:
   part=fwht([signs[x] if block>>x&1 else 0 for x in range(Q)])
   supp=[i for i,v in enumerate(part) if v];partialok &= len(supp)==q and all(abs(part[i])==q for i in supp)
   for i in supp:coverage[i]+=1
   dnow.append(tuple(sorted(i^supp[0] for i in supp)))
  union=frozenset(i for i,c in enumerate(coverage) if c)
  partialok &= all(c<=1 for c in coverage)
  if fixed is None:fixed,directions=union,dnow
  else:partialok &= fixed==union and directions==dnow
  slope=[]
  for base,tangent,_ in flats:
   for h in tangent:
    value=0
    for j,m in enumerate(genmasks):
     if p>>j&1:value^=gen(m,base^h)^gen(m,base)
    slope.append(value)
  labels.append(sum(bit<<i for i,bit in enumerate(slope)))
 rel=[z^labels[0] for z in labels];kernel={p for p,z in enumerate(rel) if not z};offs=span(blockcoords)
 assert rank(rel)==expectedslope and kernel==offs
 D=frozenset(i for i in range(Q) if not((i>>4)&1))
 assert fixed<=D and len(fixed)==expected_size and partialok
 if require_full_D:assert fixed==D
 return {"name":name,"dimension":rank(words),"phase_count":1<<len(words),"all_bent":bool(allbent),"blocks_size8":True,"nonparallel":True,"residual_size":sum(x.bit_count() for x in blocks),"fixed_total_frequency_size":len(fixed),"fixed_total_equals_u_frequency_zero_D":fixed==D,"fixed_total_nonaffine":len(fixed)&(len(fixed)-1)!=0,"partial_tiles_fixed_disjoint_pm8":bool(partialok),"translation_rank":rank(rel),"translation_kernel_offset_span":kernel==offs,"offset_rank":rank(blockcoords)}
def main():
 t=time.perf_counter()
 full=verify("four_tile_equality",[1,2,4,8,16,32,64],[A0,A1,A2,A3],[83,80,34,32],7,3,32,True)
 strict=verify("three_tile_nonaffine",[3,9,33,80],[A0,A1,A2],[9,8,5],4,1,24,False)
 out={"scope":"Exact checks of two explicit six-variable constructions only.","python":sys.version,"platform":platform.platform(),"base":"ab+cd+uv","q":q,"four_tile":full,"three_tile":strict,"elapsed_seconds":time.perf_counter()-t,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8");print(json.dumps(out,ensure_ascii=False,indent=2))
 assert full["all_bent"] and full["translation_rank"]==3 and full["fixed_total_frequency_size"]==32
 assert strict["all_bent"] and strict["translation_rank"]==1 and strict["fixed_total_frequency_size"]==24
if __name__=="__main__":main()
