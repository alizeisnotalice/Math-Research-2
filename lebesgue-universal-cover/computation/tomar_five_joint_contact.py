"""Joint rational mixtures of independent material chains at distinct anchors."""
from fractions import Fraction
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import linprog
from flint import arb
from tomar_five_contact_material import AREA,DEN,weights,q,CurvedBasin

def joint_contract(anchors,l,h):
    groups=[];nanchors=len(anchors);end=0
    for base,options in anchors:
        indexes=list(range(end,end+len(options)));end+=len(options);groups.append(indexes)
    cols=11+end+1;rows=[];rhs=[];slices=[];aslices=[];scale=2**25
    for n,(base,options) in enumerate(anchors):
        shape=[]
        for idx,event in zip(groups[n],options):
            start=len(rows)
            for v in event:
                a=np.zeros(cols);a[:11]=np.asarray(v[1:],float)*(scale/AREA);a[11+idx]=-1
                rows.append(a);rhs.append(-float(v[0])*scale/AREA)
            shape.append(slice(start,len(rows)))
        slices.append(shape)
        a=np.zeros(cols);a[:11]=np.asarray(base[1:],float)*(scale/AREA);a[11+np.asarray(groups[n],int)]=1;a[-1]=-1
        aslices.append(len(rows));rows.append(a);rhs.append(-float(base[0])*scale/AREA)
    A=np.asarray(rows);B=np.asarray(rhs);bounds=[(-1,1)]*11+[(None,None)]*(end+1)
    obj=np.zeros(cols);obj[-1]=1
    r=linprog(obj,A_ub=A,b_ub=B,bounds=bounds,method="highs")
    if not r.success:return False,l,h
    def total(dual):
        wa=weights(dual[aslices]);out=[0]*12
        for j,((base,opts),a) in enumerate(zip(anchors,wa)):
            if not a:continue
            t=[DEN*x for x in base]
            for ev,sl in zip(opts,slices[j]):
                we=weights(dual[sl])
                for k in range(12):t[k]+=sum(w*f[k]for w,f in zip(we,ev))
            for k in range(12):out[k]+=a*t[k]
        return out
    def exclude(t,L,H):
        u=sum(map(abs,t[1:]))
        if t[0]-u>0:return True
        for j,g in enumerate(t[1:]):
            if not g:continue
            d=abs(g);s=1 if g>0 else -1;N=-t[0]+u-d;x,y=int(l[j]),int(h[j])
            f=(x+y)*d+s*(y-x)*N
            if s==1:H[j]=min(int(H[j]),max(x-1,-((-f)//(2*d))))
            else:L[j]=max(int(L[j]),min(y+1,f//(2*d)))
        return bool(np.any(L>H))
    L=l.copy();H=h.copy()
    if exclude(total(-r.ineqlin.marginals),L,H):return True,L,H
    bounds[-1]=(None,0)
    for j in range(11):
        if l[j]==h[j]:continue
        for s in (-1,1):
            obj=np.zeros(cols);obj[j]=-s
            r=linprog(obj,A_ub=A,b_ub=B,bounds=bounds,method="highs")
            if not r.success or s*r.x[j]>.93:continue
            if exclude(total(-r.ineqlin.marginals),L,H):return True,L,H
    return False,L,H

class Wide(CurvedBasin):
    def __init__(self,src):
        super().__init__(src);self.scale=q(".000012")
    def contains(self,l,h):
        if np.any(np.maximum(abs(l[:2]-self.tri),abs(h[:2]-self.tri))>2**44/100000):return False
        cb=self.cb;S=2**44
        raw=[arb(arb(int(x)+int(y))/(2*S),arb(int(y)-int(x))/(2*S))for x,y in zip(l,h)]
        a=cb.unit_root(cb.add(raw[:2],cb.regular(3,0)),cb.r,1);b=cb.unit_root(cb.add(raw[:2],cb.regular(3,1)),cb.r,-1)
        if not (a[0]>0 and b[1]>0):return False
        x=((a[1]/a[0]).atan()-self.cx)/self.scale;y=(cb.pi/2-(b[0]/b[1]).atan()-self.cy)/self.scale
        if not(abs(x)<q(".5")and abs(y)<q(".5")):return False
        if self.anchors is None:
            cb.Model.degree=9;u=cb.Model({(0,0):self.cx,(1,0):self.scale});v=cb.Model({(0,0):self.cy,(0,1):self.scale})
            self.anchors=cb.model(u,v)[1]
        pose=cb.evaluate(self.anchors,x,y)
        return all(abs(raw[i]-pose[i])<q(".00007")for i in range(2,11))
