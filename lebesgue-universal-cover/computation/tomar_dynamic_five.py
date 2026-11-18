"""Membership-verified dynamic sites, integer affine fans and rational LP contraction for five bodies."""
from pathlib import Path
import numpy as np,sys,json,ctypes,math
from functools import lru_cache
from flint import arb,ctx
ctx.prec=144
from fractions import Fraction as F
P=2**40;S=2**44;goal=None;hi=None;dll=None;gen8=None
def configure(checkout,cert,library,target):
 global gen8,dll,hi,goal
 sys.path.insert(0,str(Path(checkout)/'independent'))
 import gen8 as proposal
 gen8=proposal
 exact=F(target)*2**57;goal=-(-exact.numerator//exact.denominator)
 hh=[F(z)*S for z in json.loads((Path(cert)/'part0000.json').read_bytes())['hi']]
 assert all(z.denominator==1 for z in hh)
 hi=np.array(list(map(int,hh)),dtype=np.int64)
 vec=lambda t:np.ctypeslib.ndpointer(dtype=t,ndim=1,flags='C_CONTIGUOUS')
 dll=ctypes.CDLL(str(Path(library).resolve()))
 dll.forms_for_box.argtypes=[vec(np.int64)]*4+[vec(np.int32)]*2+[ctypes.c_int,vec(np.int64)]
 dll.box_try.argtypes=[vec(np.int64)]*4+[vec(np.int32),ctypes.c_int64]
eyes=np.eye(11,dtype=np.int32);rng=np.random.default_rng(20261118);corner=rng.integers(0,2,(96,11),dtype=np.int32)*2-1; moves=np.r_[np.zeros((1,11),np.int32),eyes,-eyes,corner]
@lru_cache(50000)
def vertices(n,total):
 th=arb(total)/(2*S); R=1/(2*(arb.pi()/(2*n)).cos());v=[]
 for k in range(n):
  sn,cs=(th+2*arb.pi()*k/n).sin_cos();rx,ry=R*cs*P,R*sn*P;xx=round(float(rx.mid()));yy=round(float(ry.mid()))
  assert abs(rx-xx)<1 and abs(ry-yy)<1;v.append((xx,yy))
 return np.asarray(v,dtype=np.int64)
def table(l,h,guides=False):
 mid=(l+h)/(2*S)
 seq=gen8.hull_intervals(mid);runs=[]
 for a,b,f in seq:
  if runs and runs[-1][2]==f:runs[-1][1]=b
  else:runs.append([a,b,f])
 out=[];count=[]
 for j,n in enumerate((0,3,5,7,9)):
  idx=3*j-3;tx,ty=mid[idx:idx+2] if n else (0,0); ca=[]
  for a,b,v in runs:
   if v[5]!=j:continue
   cx,cy,r=v[2:5];cx-=tx;cy-=ty
   for u in np.linspace(a,b,1+max(1,math.ceil((b-a)/.05))) if r else (a,):
    ca.append((round((cx+r*math.cos(u))*(P-12)),round((cy+r*math.sin(u))*(P-12))))
  total=int(l[idx-1])+int(h[idx-1]) if j>=2 else 0;rot=total/(2*S)
  if guides:
   uf=32 if not n else 8*n
   ang=2*np.pi*np.arange(uf)/uf+rot
   sn=np.sin(ang);cs=np.cos(ang)
   if n:
    R=.5/np.cos(np.pi/(2*n));ids=np.rint((ang-rot+np.pi)*n/(2*np.pi)).astype(int)%n;cx=R*np.cos(rot+2*np.pi*ids/n);cy=R*np.sin(rot+2*np.pi*ids/n);pr=cs*cx+sn*cy;rad=pr+np.sqrt(1-R*R+pr*pr)
   else:rad=.5
   ca.extend(np.rint((P-12)*np.c_[cs*rad,sn*rad]).astype(np.int64).tolist())
  cand=np.unique(np.asarray(ca,np.int64).reshape(-1,2),axis=0)
  pts=cand.astype(object)
  if n:
   diff=pts[:,None,:]-vertices(n,total).astype(object);assert np.all(np.sum(diff*diff,axis=2)<(P-2)**2)
  else:assert np.all(np.sum(pts*pts,axis=1)<(P//2)**2)
  ca=np.r_[cand,np.array([[0,0],[1,0],[0,1]],np.int64)];out.append(ca.ravel());count.append(len(ca))
 return np.concatenate(out),np.asarray(count,np.int32)
def formset(l,h,m):
 allf=[]
 for dense in (False,True):
  pts,cnt=table(l,h,dense);axis=dll.box_try(l,h,hi,pts,cnt,goal)
  if axis<0:return None,axis
  fs=np.zeros((m,12),np.int64);nr=dll.forms_for_box(l,h,hi,pts,cnt,moves[:m].ravel(),m,fs.ravel());assert nr>=0;allf.extend(fs[:nr])
 return np.unique(allf,axis=0),axis

from collections import Counter
from scipy.optimize import linprog
import warnings
from scipy.optimize import OptimizeWarning
warnings.filterwarnings('ignore',category=OptimizeWarning)
OPTIONS=dict(dual_feasibility_tolerance=1e-10,primal_feasibility_tolerance=1e-10)
def combined(forms,vals):
 ids=np.flatnonzero(vals>1e-11)
 if not len(ids):return None,0
 v=vals[ids];v=v/v.max()*2**31
 ws=[max(0,round(z)) for z in v]
 return forms[ids].T.astype(object)@np.asarray(ws,dtype=object),sum(ws)
def contract(l,h,m):
 a,axis=formset(l,h,m)
 if a is None:return 1,axis,l,h,0
 if not len(a):return 0,axis,l,h,0
 g=a[:,1:].astype(float)/2**52;f=(a[:,0]-goal).astype(float)/2**52
 minmax=linprog([0]*11+[1],A_ub=np.c_[g,-np.ones(len(f))],b_ub=-f,bounds=[(-1,1)]*11+[(None,None)],method='highs',options=OPTIONS)
 if not minmax.success:return 0,axis,l,h,0
 z,den=combined(a,-minmax.ineqlin.marginals)
 if den and z[0]-sum(map(abs,z[1:]))>=den*goal:return 2,axis,l,h,0
 if minmax.fun>1e-9:return 0,axis,l,h,0
 newl=l.copy();newh=h.copy();np_cut=0
 for j in range(11):
  for sign in (1,-1):
   d=np.zeros(11);d[j]=-sign
   lp=linprog(d,A_ub=g,b_ub=-f,bounds=[(-1,1)]*11,method='highs',options=OPTIONS)
   if not lp.success or sign*lp.x[j]>.95:continue
   z,den=combined(a,-lp.ineqlin.marginals)
   if not den:continue
   Gj=z[j+1];D=sign*Gj
   if D<=0:continue
   N=goal*den-z[0]+sum(abs(val) for i,val in enumerate(z[1:]) if i!=j)
   a0=int(l[j]);b0=int(h[j]);num=(a0+b0)*D+sign*(b0-a0)*N
   if sign==1:
    cut=-((-num)//(2*D))
    if cut<b0:newh[j]=max(a0-1,int(cut));np_cut+=1
   else:
    cut=num//(2*D)
    if cut>a0:newl[j]=min(b0+1,int(cut));np_cut+=1
 if (newh[:11]<newl[:11]).any():return 3,axis,newl,newh,np_cut
 return 0,axis,newl,newh,np_cut
def work(row):
 i,label,v,limit,depth_limit=row;l=np.zeros(11,dtype=np.int64);r=l.copy();l[:11]=v[:11];r[:11]=v[11:];stack=[(l,r,0)]
 vis=0;succ=Counter();left=[];splits=0;md=0;contractions=0
 while stack and vis<limit:
  a,b,k=stack.pop();vis+=1;md=max(md,k)
  ok,axis,L,H,cut=contract(a,b,65)
  if ok:succ[ok]+=1;continue
  wid=b[:11]-a[:11];change=(wid-(H[:11]-L[:11]))*np.array([1,1,.16,1,1,.11,1,1,.09,1,1])
  eff=wid*np.array([1,1,.16,1,1,.11,1,1,.09,1,1]);ratio=1-change.max()/max(1,eff.max())
  if k>=depth_limit:left.append(np.r_[L[:11],H[:11]].tolist());continue
  if ratio<.90:
   contractions+=1;stack.append((L,H,k+1));continue
  eff=(H[:11]-L[:11])*np.array([1,1,.16,1,1,.11,1,1,.09,1,1])
  if eff[axis]*24<eff.max() or H[axis]-L[axis]<2:axis=int(np.argmax(eff))
  if H[axis]-L[axis]<2:left.append(np.r_[L[:11],H[:11]].tolist());continue
  mid=L[axis]+(H[axis]-L[axis])//2;C=L.copy();C[axis]=mid;B=H.copy();B[axis]=mid
  splits+=1;stack.extend([(C,H,k+1),(L,B,k+1)])
 left.extend(np.r_[a[:11],b[:11]].tolist() for a,b,k in stack)
 return dict(id=i,label=label,vis=vis,splits=splits,contractions=contractions,closed=dict(succ),remaining=left,md=md)
