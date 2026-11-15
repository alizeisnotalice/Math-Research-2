#!/usr/bin/env python3
"""Independent exact rectangle-union check for nonlinear cumulative/pair proof."""
from fractions import Fraction as F
from itertools import product,combinations
from functools import total_ordering
from dataclasses import dataclass
from pathlib import Path
import json,platform,hashlib
@total_ordering
@dataclass(frozen=True)
class Q:
 a:F=F(0)
 b:F=F(0)
 def __add__(s,t):
  if not isinstance(t,Q):t=Q(F(t))
  return Q(s.a+t.a,s.b+t.b)
 __radd__=__add__
 def __neg__(s):return Q(-s.a,-s.b)
 def __sub__(s,t):return s+-cv(t)
 def __rsub__(s,t):return cv(t)+-s
 def __mul__(s,t):
  t=cv(t);return Q(s.a*t.a+2*s.b*t.b,s.a*t.b+s.b*t.a)
 __rmul__=__mul__
 def sign(s):
  a,b=s.a,s.b
  if b==0:return (a>0)-(a<0)
  if a>=0 and b>0:return 1
  if a<=0 and b<0:return -1
  v=a*a-2*b*b
  return ((v>0)-(v<0)) if a>0 else ((v<0)-(v>0))
 def __lt__(s,t):return (s-cv(t)).sign()<0
 def __str__(s):return str(s.a)+'+('+str(s.b)+')sqrt2'
def cv(t):return t if isinstance(t,Q) else Q(F(t))
def union_area(rects):
 if not rects:return Q()
 xs=sorted(set([r[0] for r in rects]+[r[1] for r in rects]));area=Q()
 for left,right in zip(xs,xs[1:]):
  intervals=sorted([(r[2],r[3]) for r in rects if r[0]<=left and right<=r[1]])
  if not intervals:continue
  lo,hi=intervals[0];length=Q()
  for a,b in intervals[1:]:
   if a<=hi:hi=max(hi,b)
   else:length+=hi-lo;lo,hi=a,b
  length+=hi-lo;area+=(right-left)*length
 return area
def rectangles(bits,h,clip=True):
 R=2*h;T=8*h;out=[];sources=[(j-3)*h for j,c in enumerate(bits) if c]
 for size in range(1,len(sources)+1):
  for subset in combinations(sources,size):
   # A >=3-source subset has horizontal half-span >=2h and is excluded by open cap.
   if (max(subset)-min(subset))/2>=R:continue
   assert size<=2
   t=Q(F(3,4)*h) if size==1 else Q(F(0),F(3,4)*h)
   # Pair threshold may be too small to capture its endpoints.
   xa=cv(max(subset))-t;xb=cv(min(subset))+t
   if not xa<xb:continue
   for y in ([F(0),T] if clip else [F(0)]):
    ya=cv(y)-t;yb=cv(y)+t
    if clip:xa1=max(xa,Q());xb1=min(xb,cv(h));ya1=max(ya,Q());yb1=min(yb,cv(T))
    else:xa1,xb1,ya1,yb1=xa,xb,ya,yb
    if xa1<xb1 and ya1<yb1:out.append((xa1,xb1,ya1,yb1))
 return out

def main():
 B=Q(F(9,4),F(-3,2));assert B>Q(F(1,8))
 z=Q(F(3),F(-4,3));checks=0;zero=0
 legal=[c for c in product([0,1],repeat=7) if not any(a*b for a,b in zip(c,c[1:]))]
 vertices=[c for c in product([0,1],repeat=6) if not any(a*b for a,b in zip(c,c[1:]))]
 assert len(legal)==34 and len(vertices)==21
 for h in [F(1),F(5,7)]:
  m=F(9,4)*h*h;bh=B*h*h
  for bits in legal:
   c={j:bits[j+3] for j in range(-3,4)}
   a=union_area(rectangles(bits,h))
   formula=cv(m/2*(c[0]+c[1]))+bh*(c[-1]*c[1]+c[0]*c[2])
   assert a==formula,(h,bits,str(a),str(formula))
   delta=cv(m/2*(c[1]-c[0]))+bh*(c[-1]*c[1]-c[-2]*c[0])
   slack=delta-a+z*m*c[0]
   assert slack==bh*c[0]*(2-c[-2]-c[2]) and slack>=Q()
   if c[0]==0:assert slack==Q();zero+=1
   checks+=1
 pairbits=(0,0,1,0,1,0,0)
 pairstrip=union_area(rectangles(pairbits,F(1)))
 pairtotal=union_area(rectangles(pairbits,F(1),False))
 assert pairstrip==Q(F(27,8),F(-3,2))
 assert pairtotal==Q(F(9),F(-3))
 identities=0;caps=0
 for root in [F(0),F(1,5),F(3,4),F(1),F(2),F(5,2),F(4)]:
  m=4*root*root;t=min(F(2),root);A=4*t*t
  a={j:2*t*max(F(0),min(F(1),j+t)-max(F(0),j-t)) for j in range(-8,9)}
  assert sum(a.values(),F(0))==A
  def phi(j):return sum((v for l,v in a.items() if l>j),F(0))-(A if j<0 else F(0))
  for j in range(-6,7):
   assert phi(j-1)-phi(j)==a[j]-(A if j==0 else F(0))
   assert abs(phi(j))<=A/2<=m/2
   identities+=1
  assert phi(-3)==phi(2)==0
  if root>=2:caps+=1
 vertical_checks=0
 for T in [F(2),F(4),F(16)]:
  for m in [F(0),T*T/8,T*T/4,T*T/2]:
   for factor in [F(1),F(5,4),F(2),F(10)]:
    r=T/2*factor
    bound=m/(2*T*r)+m/(4*r*r)
    assert bound<=2*m/(T*T)<=1
    vertical_checks+=1
 for n in range(1,10):
  T=F(2**n);R=T;maxmass=F(2**n)
  assert maxmass<=T*T/2 and maxmass/4<R*R
 tail=F(9,64)+F(45,128)+F(9,64);assert tail==F(81,128)<1
 # The explicit rational rectangle uses r=33/32 and mass9/2.
 assert F(9,2)/(4*F(33,32)**2)==F(128,121)>1
 out={'python':platform.python_version(),'single_process':True,'arithmetic':'Fraction/Qsqrt2','random_seed':None,'tolerance':None,'hardcore_vertices':len(vertices),'hardcore_edges':len(legal),'complete_geometry_checks':checks,'zero_center_checks':zero,'cumulative_identities':identities,'cap_cases':caps,'single_column_periodic_tail_checks':vertical_checks,'R32_full_alphabet_maximum_checks':9,'pair_strip':str(pairstrip),'pair_total':str(pairtotal),'B':str(B),'z':str(z),'tail_at_u2':str(tail),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 dest=Path(__file__).parent/'results/nonlinear_cumulative_pair_review_exact.json';dest.parent.mkdir(exist_ok=True,parents=True);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
