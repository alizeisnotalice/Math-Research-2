#!/usr/bin/env python3
"""Corrected exact prefix search; crucially includes p_k*g_0 in each choice."""
import argparse,hashlib,json,time
from pathlib import Path
def sg(u,v):
 if u==v==0:return 0
 if u>=0 and v>=0:return 1
 if u<=0 and v<=0:return -1
 return 1 if (u*u>5*v*v if u>0 else 5*v*v>u*u) else -1
def cp(z,n):return sg(z[0]-2*n,z[1])
def ad(*z):return(sum(x for x,y in z),sum(y for x,y in z))
def iz(n):return(2*n,0)
def sc(z,n):return z[0]*n,z[1]*n
def unit(z):return cp(z,0)>=0 and cp(z,1)<=0
def fl(z):
 r=(abs(z[0])+6*abs(z[1]))//2+2;l,h=-r,r
 while l<h:
  m=(l+h+1)//2
  if cp(z,m)>=0:l=m
  else:h=m-1
 return l
def options(g,p):
 # g has entries g_0,...,g_{k-1}; include j=k term p_k*g_0.
 rest=ad(*(sc(p[j],g[-j]) for j in range(1,min(len(g),len(p)-1)+1)))
 z=fl((-rest[0],-rest[1]));return [n for n in(z,z+1) if unit(ad(iz(n),rest))]
def fcoef(g,p,k):return ad(*(sc(p[j],g[k-j]) for j in range(min(k,5)+1) if k-j<len(g)))
def km(x,y):
 u,v=x;r,w=y;return((u*r+5*v*w)//2,(u*w+v*r)//2)
def norm(p,pb):
 z=[(0,0)]*11
 for i,x in enumerate(p):
  for j,y in enumerate(pb):z[i+j]=ad(z[i+j],km(x,y))
 return tuple(u//2 for u,v in z)
def cc(g,q,h,k):return sum(h[i]*g[j]*q[k-i-j] for i in range(len(h)) for j in range(len(g)) if 0<=k-i-j<len(q))
def specs():
 for s in range(-2,3):
  for t in(-2,-1,1,2):
   for b in range(-2,3):
    for c in range(-2,3):
     for d in range(-2,3):yield s,t,b,c,d
def one(vals,D):
 s,t,b,c,d=vals;a=(2*s+t,t);ab=(2*s+t,-t);p=(iz(1),iz(b),a,iz(c),iz(d),iz(1));pb=(iz(1),iz(b),ab,iz(c),iz(d),iz(1));h=norm(p,pb);S=[((1,),(1,))];L=[]
 for k in range(1,D+1):
  N=[];pre=0
  for g,q in S:
   for x in options(g,p):
    for y in options(q,pb):
     pre+=1;gg,qq=g+(x,),q+(y,)
     if cc(gg,qq,h,k) in(0,1):N.append((gg,qq))
  L.append([k,pre,len(N)]);S=N
 free=len(S)
 for k in range(D+1,2*D+11):
  N=[];ba=bc=0
  for g,q in S:
   if not(unit(fcoef(g,p,k)) and unit(fcoef(q,pb,k))):ba+=1
   elif cc(g,q,h,k)not in(0,1):bc+=1
   else:N.append((g,q))
  L.append([k,ba,bc,len(N)]);S=N
 return {'s':s,'t':t,'b':b,'c':c,'d':d,'H':list(h),'stages':L,'free_prefix_states_after_cap':free,'post_tail_states':len(S)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('cap',type=int);ap.add_argument('--output',type=Path,required=True);ns=ap.parse_args();t=time.monotonic();rows=[one(x,ns.cap) for x in specs()]
 out={'format':'round044-x2-corrected-v1','cap':ns.cap,'processed':len(rows),'family':'s=-2..2,t=+-1,+-2,b,c,d=-2..2; P=1+b*x+a*x^2+c*x^3+d*x^4+x^5','rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-t}
 ns.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
