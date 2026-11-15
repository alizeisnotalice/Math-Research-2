#!/usr/bin/env python3
"""Independent exact full-subset strip areas in Q(sqrt2,sqrt3)."""
from fractions import Fraction as F
from dataclasses import dataclass
from functools import total_ordering
from itertools import combinations
from math import isqrt
from pathlib import Path
import json,platform,hashlib
@total_ordering
@dataclass(frozen=True)
class E:
 c:tuple=(F(0),F(0),F(0),F(0))
 def __add__(s,t):
  t=cv(t);return E(tuple(a+b for a,b in zip(s.c,t.c)))
 __radd__=__add__
 def __neg__(s):return E(tuple(-a for a in s.c))
 def __sub__(s,t):return s+-cv(t)
 def __rsub__(s,t):return cv(t)+-s
 def __mul__(s,t):
  t=cv(t);out=[F(0)]*4
  for i,a in enumerate(s.c):
   for j,b in enumerate(t.c):
    common=i&j;factor=(2 if common&1 else 1)*(3 if common&2 else 1)
    out[i^j]+=a*b*factor
  return E(tuple(out))
 __rmul__=__mul__
 def sign(s):
  if all(a==0 for a in s.c):return 0
  for bits in [8,16,32,64,128,256,512,1024]:
   D=2**bits;lo=hi=s.c[0]
   for coef,rad in zip(s.c[1:],[2,3,6]):
    floor=isqrt(rad*D*D);a=F(floor,D);b=F(floor+1,D)
    lo+=coef*(a if coef>=0 else b);hi+=coef*(b if coef>=0 else a)
   if lo>0:return 1
   if hi<0:return -1
  raise ArithmeticError('uncertified sign')
 def __lt__(s,t):return (s-cv(t)).sign()<0
 def __str__(s):return '['+','.join(str(a) for a in s.c)+'] basis(1,sqrt2,sqrt3,sqrt6)'
def cv(t):return t if isinstance(t,E) else E((F(t),F(0),F(0),F(0)))
def alg(a=0,b=0,c=0,d=0):return E(tuple(F(v) for v in [a,b,c,d]))
def area(rects):
 if not rects:return alg()
 xs=sorted(set([r[0] for r in rects]+[r[1] for r in rects]));total=alg()
 for left,right in zip(xs,xs[1:]):
  ys=sorted([(r[2],r[3]) for r in rects if r[0]<=left and right<=r[1]])
  if not ys:continue
  lo,hi=ys[0];length=alg()
  for a,b in ys[1:]:
   if a<=hi:hi=max(hi,b)
   else:length+=hi-lo;lo,hi=a,b
  total+=(right-left)*(length+hi-lo)
 return total
points=[(F(-1),F(0)),(F(1),F(3,2)),(F(2),F(0))]
roots={1:alg(1),2:alg(0,1),3:alg(0,0,1)}
def rectangles(config):
 out=[];raw=0
 for count in range(1,len(config)+1):
  for sub in combinations(config,count):
   raw+=1;t=roots[count]
   xa=cv(max(points[i][0] for i in sub))-t;xb=cv(min(points[i][0] for i in sub))+t
   ya=cv(max(points[i][1] for i in sub))-t;yb=cv(min(points[i][1] for i in sub))+t
   # Different period-label atoms cannot mix for r<2 (separation13/2>4).
   for lift in [-1,0,1]:
    a=max(xa,cv(0));b=min(xb,cv(F(1,2)));c=max(ya+8*lift,cv(0));d=min(yb+8*lift,cv(8))
    if a<b and c<d:out.append((a,b,c,d))
 return out,raw

def main():
 areas={};all_rects={};raw=0
 for size in range(1,4):
  for sub in combinations(range(3),size):
   rects,n=rectangles(sub);raw+=n;all_rects[sub]=rects;areas[sub]=area(rects)
 assert raw==19 and len(areas)==7
 assert areas[(0,)]==alg() and areas[(1,)]==alg(1) and areas[(2,)]==alg()
 assert areas[(0,1)]==alg(4,-2) and areas[(0,2)]==alg() and areas[(1,2)]==alg(1)
 pred=sum((areas[s] for s in combinations(range(3),2)),alg())-sum((areas[(i,)] for i in range(3)),alg())
 full=areas[(0,1,2)];res=full-pred
 assert full==alg(F(7,2),2,F(-3,2),-1)
 assert res==alg(F(-1,2),4,F(-3,2),-1)
 pure=area([r for sub,rects in all_rects.items() if len(sub)<=2 for r in rects])
 assert pure==pred==alg(4,-2)
 assert res>cv(F(211,1980))>cv(F(1,10))
 assert F(140,99)**2<2 and F(26,15)**2>3 and F(49,20)**2>6
 assert -F(1,2)+4*F(140,99)-F(3,2)*F(26,15)-F(49,20)==F(211,1980)
 # Explicit U closed-corner extrema retain the strict three-source witness.
 corners=[(x,y) for x in [F(7,16),F(1,2)] for y in [F(-1,32),F(1,32)]]
 assert all(max(abs(x-a),abs(y-b))<F(13,8) for x,y in corners for a,b in points)
 assert 12/(4*F(13,8)**2)==F(192,169)>1
 assert F(23,16)**2>2 and F(3,2)**2>2
 assert 12/(4*F(2)**2)==F(3,4)<1
 assert F(3)/(4*F(13,4))+F(3)/F(13,4)**2==F(87,169)<1
 # Arithmetic behind the range lower bound, exact quarter-integer z cases.
 rangechecks=0
 for L in [1,2,3,5]:
  for n in range(1,25):
   p=-1
   while 16*L*L*4**(p+1)<=2**n:p+=1
   if p>=0:
    z=F(p+1,4)
    assert 16*L*L*2**int(8*z)>2**n
   else:
    z=F(0);assert 16*L*L>2**n
   D=F(L)-F(1,2)
   assert L>D and L<=D+1
   assert 16*(D+1)**2*2**int(8*z)>2**n
   rangechecks+=1
 out={'python':platform.python_version(),'single_process':True,'arithmetic':'Fraction/exact Qsqrt2sqrt3, sign by directed isqrt enclosure','random_seed':None,'tolerance':None,'source_configurations':len(areas),'raw_nonempty_subsets':raw,'areas':{str(k):str(v) for k,v in areas.items()},'pair_prediction':str(pred),'purepair_union':str(pure),'residual':str(res),'strict_rational_residual_lower_bound':str(F(211,1980)),'range_growth_arithmetic_checks':rangechecks,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 dest=Path(__file__).parent/'results/pair_cumulative_triple_review_exact.json';dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
