import sys,json
import numpy as np
from scipy.optimize import linprog
from pathlib import Path
N=int(sys.argv[1]) if len(sys.argv)>1 else 31
ks=np.array([k for k in range(5,N+1,2) if k%3]);signs=np.where(ks%6==1,1,-1)
def mat(t,curve=False):
 t=np.array(t);T=t[:,None]*ks
 if curve:
  return np.column_stack([np.sin(3*t),np.cos(T),np.sin(T)])
 return np.column_stack([-np.sin(3*t)/8,
     (np.cos(T)-np.cos(t)[:,None])/(1-ks**2),
     (np.sin(T)-np.sin(t)[:,None]*signs)/(1-ks**2)])
rads=np.pi/180
spec=[(31.50,[1,1,-1]),(31.54,[1,-1,-1]),(28.50,[1,1,-1]),(28.46,[-1,1,-1])]
if len(sys.argv)>4 and sys.argv[4]=='reverse':
 spec=[(31.54,[1,1,-1]),(31.50,[1,-1,-1]),(28.46,[1,1,-1]),(28.50,[-1,1,-1])]
roll=next((int(a.split('=')[1]) for a in sys.argv[4:] if a.startswith('roll=')),0)
if roll:
 pa,pb=np.roll(spec[0][1],roll),np.roll(spec[1][1],roll)
 spec=[(spec[0][0],pa),(spec[1][0],pb),(spec[2][0],pa[[1,0,2]]),(spec[3][0],pb[[1,0,2]])]
theta=[(t+120*j)*rads for t,l in spec for j in range(3)];sgn=np.array([c for t,l in spec for c in l])
A=mat(theta)*sgn[:,None]
phi=(np.arange(1440)+.37)*2*np.pi/1440
C=mat(phi,True)
focus=int(sys.argv[3]) if len(sys.argv)>3 else 0
coef=np.ones(12) if not focus else np.array([int((3 if j//3 in (0,2) else 5)==focus) for j in range(12)])
con=np.vstack([np.column_stack([-A,coef]),np.column_stack([C,np.zeros(len(C))]),np.column_stack([-C,np.zeros(len(C))])])
d=np.r_[-(1-coef)*1e-8,np.full(2*len(C),.499)]
obj=np.zeros(con.shape[1]);obj[-1]=-1
res=linprog(obj,A_ub=con,b_ub=d,bounds=([(.000001,None)] if 'negroot' in sys.argv[4:] else [(None,-.000001)])+[(None,None)]*(len(obj)-1),method='highs')
print(res.success,res.message)
if res.success:
 x=res.x;print('margin',x[-1],'ratio',x[-1]/np.deg2rad(.04),'H3',x[0], 'curv Lipschitz',np.sum(np.abs(x[1:1+len(ks)])*ks)+np.sum(np.abs(x[1+len(ks):-1])*ks)+3*abs(x[0]))
 data=dict(focus=focus,max_degree=N,frequencies=ks.tolist(),coefficients=[f'{v:.13f}' for v in x[:-1]],min_common_support=float(x[-1]))
 Path(sys.argv[2]).write_text(json.dumps(data,indent=1)+'\n')
