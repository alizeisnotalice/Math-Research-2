"""Test contact-chart cells with all 32768 available chains factored by their six independent events."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json,sys
from collections import Counter
from pathlib import Path
import numpy as np
from flint import arb,ctx
from scipy.optimize import linprog
from scipy import sparse

class Chart:
    def __init__(self,source,cx,cy,radius,degree,d,cube):
        sys.path.insert(0,str(source.resolve()))
        import certify_five_body_contact_basin as b
        import certify_five_body_contact_transverse as tr
        self.b=b;self.tr=tr;ctx.prec=192
        self.d=b.q(d);self.rho=Fraction(cube);self.scale=b.q(radius)
        b.Model.degree=degree
        x=b.Model({(0,0):b.q(cx),(1,0):self.scale})
        y=b.Model({(0,0):b.q(cy),(0,1):self.scale})
        self.data=b.model(x,y);self.vol=b.area(x,y)
        self.sign=np.array(list(product((-1,1),repeat=9)),np.int64)
    def grad_floor(self,base,all_ev):
        # A single affine probability rule proposes all 512 sign combinations.
        # The rounded weights below are independently checked as probabilities.
        widths=(8,8,8,8,4,2);n=sum(widths);M=np.zeros((6,n));offsets=[0]
        for j,w in enumerate(widths):
            M[j,offsets[-1]:offsets[-1]+w]=1;offsets.append(offsets[-1]+w)
        G=np.array([[float(a.mid())for a in row]for row in all_ev]);b=np.array([float(a.mid())for a in base])
        H=sparse.csc_array(np.r_[M,G.T]);E=sparse.block_diag([H]*10,format="csc")
        D=sparse.csc_array((150,9*n))
        rad=np.zeros((10,15));rad[1:,6:]=-np.eye(9)
        Aeq=sparse.hstack([E,D,sparse.csc_array(rad.ravel()[:,None])],format="csc")
        beq=np.zeros((10,15));beq[0,:6]=1;beq[0,6:]=-b
        I=sparse.eye_array(9*n,format="csc");Z=sparse.csc_array((9*n,n));I0=sparse.eye_array(n,format="csc")
        up=sparse.vstack([sparse.hstack([Z,I,-I]),sparse.hstack([Z,-I,-I]),
            sparse.hstack([-I0,sparse.csc_array((n,9*n)),sparse.hstack([I0]*9)])])
        up=sparse.hstack([up,sparse.csc_array((up.shape[0],1))],format="csc")
        obj=np.r_[np.zeros(up.shape[1]-1),-1.]
        result=linprog(obj,A_eq=Aeq,b_eq=beq.ravel(),A_ub=up,b_ub=np.zeros(up.shape[0]),bounds=[(None,None)]*len(obj),method="highs")
        assert result.success
        guess=result.x[:n]+self.sign@result.x[n:10*n].reshape(9,n)
        W=np.zeros((512,n),dtype=np.int64)
        for i in range(6):
            ix=slice(offsets[i],offsets[i+1]);vv=np.maximum(0,guess[:,ix]);assert np.all(vv.sum(axis=1)>0)
            candidate=np.floor(vv/vv.sum(axis=1,keepdims=True)*2**34).astype(np.int64)
            candidate[np.arange(512),vv.argmax(axis=1)]+=2**34-candidate.sum(axis=1)
            assert np.all(candidate>=0) and np.all(candidate.sum(axis=1)==2**34)
            W[:,ix]=candidate
        grid=2**24;den=2**34;v=[base,*all_ev]
        # Each of six groups is a probability vector of mass den. With all absolute
        # entries under one, every integer dot product here has magnitude below 7*2^58.
        assert all(abs(y)<1 for row in v for y in row)
        lo=np.array([[int((y.lower()*grid).floor().unique_fmpz())for y in row]for row in v],np.int64)
        hi=np.array([[int((y.upper()*grid).ceil().unique_fmpz())for y in row]for row in v],np.int64)
        inf=lo[0]*den+W@lo[1:];sup=hi[0]*den+W@hi[1:]
        margin=int(np.where(self.sign>0,inf,-sup).min())
        return Fraction(margin,den*grid)
    def test(self,x,y):
        b=self.b;tr=self.tr;d=self.d
        raw=b.geometry_box(self.data,x,y);pcs,t,points,rots=raw;b.atomic_order(pcs,rots)
        H=[[self.vol.second(i,j,x,y,self.scale)-b.q(".1")*(i==j)for j in range(2)]for i in range(2)]
        assert H[0][0]>0 and b.det(*H)>0,("hessian",str(b.det(*H)))
        t,links,arcs,baseline=b.links_and_arcs(raw)
        si,co=arb(0,d).sin_cos()
        def matrix():return [[arb(0)for _ in range(9)]for _ in range(9)]
        def addmat(A,B):return [[v+w for v,w in zip(a,b)]for a,b in zip(A,B)]
        def moved(i,P):
            if i in (0,3):return P,{},{}
            k=3*((i-5)//2);pos=t[k+3:k+5];v=b.rotated(b.sub(P,pos),co,si)
            return b.add(b.add(pos,(arb(0,d),arb(0,d))),v),{k:b.left(v),k+1:(1,0),k+2:(0,1)},{k:b.scale(-1,v)}
        HH=matrix();ori=[]
        for owner,span in arcs:
            k=3*((owner-5)//2);pose=b.add(t[k+3:k+5],(arb(0,d),arb(0,d)));v=b.rotated(span,co,si)
            HH[k][k]-=b.det(pose,v)/2
            for j,e in enumerate(((1,0),(0,1)),k+1):HH[k][j]+=b.det(e,b.left(v))/2;HH[j][k]=HH[k][j]
        def edge(a,c):
            mat=matrix();i,P=a;j,Q=c;V,dv,hv=moved(i,P);W,dw,hw=moved(j,Q)
            ori.append(b.det(V,W).lower())
            for k,dp in hv.items():mat[k][k]+=b.det(dp,W)/2
            for k,dq in hw.items():mat[k][k]+=b.det(V,dq)/2
            for k,dp in dv.items():
                for j,dq in dw.items():
                    z=b.det(dp,dq)/2;mat[k][j]+=z;mat[j][k]+=z
            return tr.connector(*a,*c,t),mat
        groups={}
        for tag,seq in links:
            if tag not in tr.EVENT_ORDER:
                g,h=edge(*seq);HH=addmat(HH,h);baseline=tr.plus(baseline,g);continue
            alts=[];ms=[]
            for z in product((False,True),repeat=len(seq)-2):
                v=[seq[0],*(u for u,s in zip(seq[1:-1],z)if s),seq[-1]]
                grad=[arb(0)]*9;mat=matrix()
                for a,c in zip(v,v[1:]):
                    gg,hh=edge(a,c);grad=tr.plus(grad,gg);mat=addmat(mat,hh)
                alts.append(grad);ms.append(mat)
            union=ms[0]
            for h in ms[1:]:union=[[u.union(v)for u,v in zip(a,c)]for a,c in zip(union,h)]
            HH=addmat(HH,union);groups[tag]=alts
        margin=self.grad_floor(baseline,[v for tag in tr.EVENT_ORDER for v in groups[tag]])
        assert min(ori)>0 and margin>self.rho,("gradient",str(margin),str(min(ori)))
        for i,n in enumerate((3,5,7,9)):
            assert abs(t[3*i])+abs(t[3*i+1])+2*d<1-1/(2*(b.pi/(2*n)).cos())
        tail=max(sum(abs(z).abs_upper()for z in row)for row in HH)*d/2
        lin=b.q(margin)-tail
        assert lin>0,("tail",str(margin),str(tail))
        return b.det(*H).lower(),margin,min(ori),lin.lower()

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source",type=Path);p.add_argument("--adaptive-max",type=int)
    p.add_argument("--radius",default=".000024");p.add_argument("--transverse",default=".0001");p.add_argument("--cube",default=".00014")
    p.add_argument("--degree",type=int,default=9);p.add_argument("--cells",type=int,default=30);p.add_argument("--only")
    a=p.parse_args();ctx.prec=192;center=(".513440833040237835","1.578292413856733933")
    c=Chart(a.source,*center,a.radius,a.degree,a.transverse,a.cube)
    initial=a.cells;results=[];counts=Counter();nodes=[]
    z=range(initial**2)if not a.only else list(map(int,a.only.split(",")))
    for k in z:
        i,j=divmod(k,initial);assert 0<=i<initial and 0<=j<initial
        nodes.append((initial,i,j))
    while nodes:
        n,i,j=nodes.pop()
        unit=lambda v:arb(c.b.q(Fraction(2*v+1-n,2*n)),c.b.q(Fraction(1,2*n)))
        try:m=c.test(unit(i),unit(j))
        except AssertionError as err:
            if a.adaptive_max is None or 2*n>a.adaptive_max:raise AssertionError((n,i,j,err.args))from err
            nodes.extend((2*n,2*i+u,2*j+v)for u in (0,1)for v in (0,1))
            continue
        results.append(m);counts[n]+=1
    assert sum(Fraction(num,n*n) for n,num in counts.items())==Fraction(len(z),initial**2)
    print(json.dumps(dict(center=center,scale=a.radius,halfwidth=str(Fraction(a.radius)/2),
        transverse=a.transverse,cube=a.cube,degree=a.degree,initial=initial,accepted_cells=counts,total_checked=len(results),complete=not a.only,
        minima=[str(min(x[k]for x in results))for k in range(4)],error=str(c.vol.error),
        code_sha256=sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))

if __name__=="__main__":main()
