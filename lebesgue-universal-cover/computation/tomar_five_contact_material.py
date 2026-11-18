"""Affine support of the exact material chains at a two-angle five-body anchor."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import sys
import numpy as np
from flint import arb,ctx
from scipy.optimize import linprog

S=2**44
AREA=2**68
DEN=2**32
def rot(v,s,c):return (c*v[0]-s*v[1],s*v[0]+c*v[1])
def add(p,q):return (p[0]+q[0],p[1]+q[1])
def sub(p,q):return (p[0]-q[0],p[1]-q[1])
def sc(s,p):return (s*p[0],s*p[1])
def J(p):return (-p[1],p[0])
def det(p,q):return p[0]*q[1]-p[1]*q[0]
def norm(p):return abs(p[0])+abs(p[1])
def zz():return (arb(0),arb(0))
def weights(xs):
    v=[max(0,float(x))for x in xs];total=sum(v)
    if not total:return [DEN]+[0]*(len(v)-1)
    w=[int(x/total*DEN)for x in v];w[v.index(max(v))]+=DEN-sum(w)
    assert min(w)>=0 and sum(w)==DEN
    return w
def quant(row):
    v=row[0];out=[]
    for t in row[1:]:
        a=round(float(t.mid())*AREA)
        v-=abs(t-arb(a)/AREA);out.append(a)
    h=int((v.lower()*AREA).floor().unique_fmpz())
    return [h,*out]
def q(fr):
    v=Fraction(fr)
    return arb(v.numerator)/v.denominator

class Material:
    def __init__(self,source):
        sys.path.insert(0,str(Path(source).resolve()))
        import certify_five_body_contact_stationary as st
        import certify_five_body_contact_transverse as tr
        self.st=st;self.tr=tr;ctx.prec=192
        self.fixed=tuple(st.interval(*v)for v in st.BOX)
        self.fixed_anchor=self.anchor(*self.fixed)
    def anchor(self,a,b):
        st=self.st;tr=self.tr
        st.verify_topology(a,b)
        pcs,t,p,_=st.model(a,b)
        us={tag:st.normal(end)for tag,end,*_ in pcs if tag in tr.EVENT_ORDER}
        order=tr.face_orders(p,us);start=arb(0);arcs=[];fixed=[];events={}
        for p1,p2 in zip(pcs,pcs[1:]+pcs[:1]):
            tag,end,i,_,rad,v=p1;_,_,j,_,rad2,w=p2
            n0,n1=st.normal(start),st.normal(end)
            if rad and i:arcs.append((i,sc(rad,sub(n1,n0))))
            P=(i,add(v,sc(rad,n1)));Q=(j,add(w,sc(rad2,n1)))
            if i!=j:
                if tag in order:
                    seq=order[tag];assert seq[0][0]==i and seq[-1][0]==j
                    chain=[]
                    for z in product((False,True),repeat=len(seq)-2):
                        vv=[seq[0],*(u for u,y in zip(seq[1:-1],z)if y),seq[-1]]
                        chain.append(list(zip(vv,vv[1:])))
                    events[tag]=chain
                else:fixed.append((P,Q))
            start=end
        return t,st.area(a,b),arcs,fixed,[events[k]for k in tr.EVENT_ORDER]
    def near_anchor(self,l,h):
        st=self.st
        xy=tuple(arb(int(l[i])+int(h[i]))/(2*S)for i in (0,1))
        A=st.unit_root(add(xy,st.regular(3,0)),st.r,1)
        B=st.unit_root(add(xy,st.regular(3,1)),st.r,-1)
        assert A[0]>0 and B[1]>0
        return self.anchor((A[1]/A[0]).atan(),st.pi/2-(B[0]/B[1]).atan())
    def formsets(self,l,h,anchor,target):
        pose,M,arcs,fixed,events=anchor
        center=[arb(int(a)+int(b))/(2*S)for a,b in zip(l,h)]
        widths=[arb(int(b)-int(a))/(2*S)for a,b in zip(l,h)]
        motion={0:(zz(),zz(),arb(0),arb(1),None,arb(0),arb(0),arb(0),arb(0))}
        for body,k in ((3,0),(5,3),(7,6),(9,9)):
            now=tuple(center[k:k+2]);old=tuple(pose[k:k+2]);idx=k-1 if k else None
            da=center[idx]-pose[idx] if k else arb(0);wr=widths[idx]if k else arb(0)
            s,c=da.sin_cos();dt=widths[k]+widths[k+1]
            motion[body]=(now,old,s,c,idx,wr,dt,abs(da)+wr,norm(sub(now,old))+dt)
            if not norm(now)+dt<1-1/(2*(self.st.pi/(2*body)).cos()):return None
        def point(owner,p):
            if owner==0:return p,{},arb(0),arb(0),arb(0)
            now,old,s,c,idx,wr,dt,fullrot,fulltr=motion[owner]
            v=rot(sub(p,old),s,c);vsize=norm(sub(p,old));k=0 if owner==3 else idx+1
            deriv={k:(arb(1),arb(0)),k+1:(arb(0),arb(1))}
            if idx is not None:deriv[idx]=J(v)
            return add(now,v),deriv,dt+vsize*wr,vsize*wr*wr/2,fulltr+vsize*fullrot
        def chord(p,q):
            a,du,dp,ep,da=point(*p);b,dv,dq,eq,db=point(*q)
            ori=det(p[1],q[1])-norm(p[1])*db-norm(q[1])*da-da*db
            if not ori>0:return None
            error=(norm(a)*eq+norm(b)*ep+dp*dq)/2
            value=(det(a,b)-det(p[1],q[1]))/2-error
            return [value,*[(det(du.get(i,zz()),b)+det(a,dv.get(i,zz())))*widths[i]/2 for i in range(11)]]
        def plus(v,w):return [x+y for x,y in zip(v,w,strict=True)]
        base=[M-q(target),*[arb(0) for _ in range(11)]]
        for owner,span in arcs:
            t,origin,s,c,idx,wr,dt,_,_=motion[owner];k=0 if owner==3 else idx+1
            v=rot(span,s,c);error=norm(span)*(dt*wr+norm(t)*wr*wr/2)/2
            base[0]+=(det(t,v)-det(origin,span))/2-error
            base[1+k]+=v[1]*widths[k]/2;base[2+k]-=v[0]*widths[k+1]/2
            if idx is not None:base[1+idx]+=det(t,J(v))*wr/2
        for p,q1 in fixed:
            v=chord(p,q1)
            if v is None:return None
            base=plus(base,v)
        terms=[]
        for event in events:
            opts=[]
            for path in event:
                value=[arb(0)]*12
                for p,q1 in path:
                    v=chord(p,q1)
                    if v is None:break
                    value=plus(value,v)
                else:opts.append(quant(value))
            if not opts:return None
            terms.append(opts)
        return quant(base),terms

def solve(forms):
    base,events=forms
    widths=list(map(len,events));n=len(events);cols=11+n
    obj=np.r_[np.asarray(base[1:],float)/AREA,np.ones(n)]
    mats=[];rhs=[]
    for j,ev in enumerate(events):
        for v in ev:
            a=np.zeros(cols);a[:11]=np.asarray(v[1:],float)/AREA;a[11+j]=-1;mats.append(a);rhs.append(-v[0]/AREA)
    # Scale the area units to make very small gap constraints visible to HiGHS.
    scale=2**25;A=np.asarray(mats);A[:,:11]*=scale
    obj[:11]*=scale
    res=linprog(obj,A_ub=A,b_ub=np.asarray(rhs)*scale,bounds=[(-1,1)]*11+[(None,None)]*n,method="highs")
    if not res.success:return None
    dual=-res.ineqlin.marginals;total=[DEN*x for x in base];off=0
    for ev,w in zip(events,widths):
        ws=weights(dual[off:off+w]);off+=w
        for i in range(12):total[i]+=sum(k*v[i]for k,v in zip(ws,ev))
    margin=total[0]-sum(map(abs,total[1:]))
    return margin,total

def contract_forms(forms,l,h):
    """Use LP only to propose exact per-event convex weights and outward cuts."""
    base,events=forms
    n=len(events);m=list(map(len,events));A=[];B=[]
    scale=2**25
    for k,ev in enumerate(events):
        for v in ev:
            row=[float(g)*scale/AREA for g in v[1:]]+[0.]*n
            row[11+k]=-1;A.append(row);B.append(-float(v[0])*scale/AREA)
    c=np.r_[np.asarray(base[1:],float)*scale/AREA,np.ones(n)]
    A=np.asarray(A);B=np.asarray(B);bounds=[(-1,1)]*11+[(None,None)]*n
    first=linprog(c,A_ub=A,b_ub=B,bounds=bounds,method="highs")
    if not first.success:return False,l,h
    def pooled(dual):
        total=[DEN*x for x in base];off=0
        for ev,num in zip(events,m):
            z=weights(dual[off:off+num]);off+=num
            for i in range(12):total[i]+=sum(w*v[i]for w,v in zip(z,ev))
        return total
    def shrink(t,L,H):
        if t[0]-sum(map(abs,t[1:]))>0:return True
        normg=sum(map(abs,t[1:]))
        for j,g in enumerate(t[1:]):
            if g==0:continue
            den=abs(g);sign=1 if g>0 else -1
            N=-t[0]+normg-den
            lo,hi=int(l[j]),int(h[j])
            num=(lo+hi)*den+sign*(hi-lo)*N
            if sign==1:H[j]=min(int(H[j]),max(lo-1,-((-num)//(2*den))))
            else:L[j]=max(int(L[j]),min(hi+1,num//(2*den)))
        return bool(np.any(L>H))
    L=l.copy();H=h.copy()
    if shrink(pooled(-first.ineqlin.marginals),L,H):return True,L,H
    A=np.vstack([A,c]);B=np.r_[B,-float(base[0])*scale/AREA]
    for j in range(11):
        if h[j]==l[j]:continue
        for sign in (-1,1):
            obj=np.zeros(11+n);obj[j]=-sign
            result=linprog(obj,A_ub=A,b_ub=B,bounds=bounds,method="highs")
            if not result.success or sign*result.x[j]>.95:continue
            t=pooled(-result.ineqlin.marginals[:-1])
            if shrink(t,L,H):return True,L,H
    return False,L,H

class CurvedBasin:
    def __init__(self,source):
        sys.path.insert(0,str(Path(source).resolve()))
        import certify_five_body_contact_basin as cb
        self.cb=cb;self.anchors=None
        self.cx,self.cy=map(q,cb.CENTERS);self.scale=q("0.000006")
        self.tri=np.asarray([round(float(q(v))*S)for v in cb.RAW[:2]],np.int64)
    def contains(self,l,h):
        # This cheap guard selects candidates for the exact curved-basin test.
        if np.any(np.maximum(abs(l[:2]-self.tri),abs(h[:2]-self.tri))>S/100000):return False
        cb=self.cb
        raw=[arb(arb(int(x)+int(y))/(2*S),arb(int(y)-int(x))/(2*S))for x,y in zip(l,h)]
        a=cb.unit_root(cb.add(raw[:2],cb.regular(3,0)),cb.r,1)
        b=cb.unit_root(cb.add(raw[:2],cb.regular(3,1)),cb.r,-1)
        if not (a[0]>0 and b[1]>0):return False
        x=((a[1]/a[0]).atan()-self.cx)/self.scale
        y=(cb.pi/2-(b[0]/b[1]).atan()-self.cy)/self.scale
        if not (abs(x)<q(".5") and abs(y)<q(".5")):return False
        if self.anchors is None:
            cb.Model.degree=9
            u=cb.Model({(0,0):self.cx,(1,0):self.scale})
            v=cb.Model({(0,0):self.cy,(0,1):self.scale})
            self.anchors=cb.model(u,v)[1]
        pose=cb.evaluate(self.anchors,x,y)
        return all(abs(raw[i]-pose[i])<q(".00008") for i in range(2,11))
