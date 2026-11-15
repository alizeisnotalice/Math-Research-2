#!/usr/bin/env python3
"""Reconstruct from finite W; no generator imports; exact tensor Bernstein."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json,sys,hashlib
import sympy as S
s,u=S.symbols('s u'); h=1+s

def build(k,R):
 W=lambda j:sum((i+1)*(-s)**i for i in range(j+1))
 w=W(k); un=u/h; unm=(s*s-1+u)/(h*s*s); unp=(1-s*s+s*s*u)/h
 fm=S.cancel(h**3*s*s*(w*unm-un*W(k-2)))
 fp=S.cancel(h**3*(w*unp-un*W(k+2)))
 z=1-u
 q=(2*s*s-s)*h*h*w*w+(1-z)*(1-s-2*s*s*z)*w+(1-z)**2*(s*W(k+1)-W(k+2))
 return {'g1':s**(R+k+1)-1+u,'g2':s**(R+k+1)+1-u,'g3h':h*h*w-h*u,'fm_h2':fm,'fp_h2':fp,'F_h2q':S.expand(h*h*q)}

def bounds(expr,box):
 p=S.Poly(expr,s,u); nx,ny=p.degree(s),p.degree(u); a,b,c,d=box
 power={}
 for (j,k),v in p.terms():
  v=Q(int(v.p),int(v.q))
  for i in range(j+1):
   for l in range(k+1):
    power[i,l]=power.get((i,l),Q(0))+v*comb(j,i)*a**(j-i)*(b-a)**i*comb(k,l)*c**(k-l)*(d-c)**l
 bs=[sum(v*Q(comb(i,j),comb(nx,j))*Q(comb(l,k),comb(ny,k)) for (j,k),v in power.items() if j<=i and k<=l) for i in range(nx+1) for l in range(ny+1)]
 return min(bs),max(bs)

def verify(path):
 rec=json.loads(path.read_text()); k,R=rec['k'],rec['r']; assert k in (4,6) and R==k+3
 root=(Q(1,4*k*(k+1)),Q(1),Q(0),Q(2)); assert tuple(map(Q,rec['root']))==root
 polys=build(k,R); assert S.expand(polys['fm_h2'].subs(s,1))==4*u; assert S.expand(polys['fp_h2'].subs(s,1))==-4*u
 boxes=[]
 for item in rec['leaves']:
  box=tuple(map(Q,item['box'])); a,b,c,d=box
  assert root[0]<=a<b<=root[1] and root[2]<=c<d<=root[3]
  leaf=item['leaf']; name=leaf.get('guard','F_h2q'); lo,hi=bounds(polys[name],box)
  if leaf['kind']=='guard_upper_negative': assert name!='F_h2q' and hi<0 and hi==Q(leaf['bound'])
  else: assert leaf['kind']=='F_lower_positive' and lo>0 and lo==Q(leaf['bound'])
  boxes.append(box)
 xs=sorted({root[0],root[1],*(v for b in boxes for v in b[:2])}); ys=sorted({root[2],root[3],*(v for b in boxes for v in b[2:])})
 cells=0
 for a,b in zip(xs,xs[1:]):
  for c,d in zip(ys,ys[1:]):
   assert any(x<=a and b<=y and z<=c and d<=w for x,y,z,w in boxes)
   cells+=1
 return {'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'k':k,'R':R,'leaves':len(boxes),'closed_cells_covered':cells,'status':'verified_from_finite_W'}
if __name__=='__main__':
 print(json.dumps([verify(Path(p)) for p in sys.argv[1:]],indent=2))
