"""Certify the signed split of a third strip by bounded equilibria and two harmonics."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from math import cos, sin, pi
from pathlib import Path
import argparse,json,sys
import numpy as np
from flint import arb,ctx,__version__ as flint_version
ctx.prec=160

def q(z):
    a=F(z);return arb(a.numerator)/a.denominator
def deg(z):return q(z)*arb.pi()/180
def unit(x):return x.cos(),x.sin()
def cross(x,y):return x[0]*y[1]-x[1]*y[0]
def vsum(coords,w):
    return tuple(sum(v[j]*x for v,x in zip(coords,w)) for j in (0,1))
def enc(x):
    m,k=x.lower().man_exp()
    return str(F(int(m))*F(2)**int(k))

def directions(d):
    return [deg(0),deg(64),deg(96),deg(128)-d,deg(180),deg(212),deg(244),deg(308)+d]
def gaps(T):
    return [v-u for u,v in zip(T,T[1:]+[T[0]+2*arb.pi()])]
def patterns(t,L):
    normals=np.c_[np.cos(t),np.sin(t)];ell=np.array(L)
    for i,j in combinations(range(8),2):
        f=[k for k in range(8) if k not in (i,j)]
        mat=normals[[i,j]].T
        if abs(np.linalg.det(mat))<.15:continue
        masks=np.array(list(product((0,1),repeat=6)))
        guess=ell[f]*masks
        v=np.linalg.solve(mat,(-guess@normals[f]).T).T
        for mask in masks[np.all((v>=-1e-6)&(v<=ell[[i,j]]+1e-6),axis=1)]:
            yield i,j,f,mask.tolist()

def run(lo,hi,er,n,free,target):
    target=q(target);B=q('.844082046521218581802948589770767704354')
    phase=deg(144);blend=q('.01');s=q('.999');root3=arb(3).sqrt()
    results=[]
    for panel in range(n):
        a=lo+(hi-lo)*panel/n;b=lo+(hi-lo)*(panel+1)/n
        center=deg((a+b)/2);rad=deg((b-a)/2);e=deg(er)
        mids=directions(center);box=deg(a).union(deg(b))
        T=[v+(-e).union(e) for v in directions(box)];g=gaps(T)
        N=list(map(unit,T));L=[((g[i-1]/2).tan()+(g[i]/2).tan())/2 for i in range(8)]
        assert min(L)>0 and min(g)>0 and max(g)<arb.pi()/2
        partial=[((g[i-1]/2).cos()**(-2)-(g[i]/2).cos()**(-2))/8 for i in range(8)]
        A=sum((x/2).tan() for x in gaps(mids))/4+(partial[7]-partial[3])*(-rad).union(rad)+sum(z*(-e).union(e) for z in partial)
        gamma=max(((1/g[i-1].tan()+1/g[i].tan())/L[i]**2).upper() for i in range(8))
        if A>target:
            results.append(dict(panel=panel,d=[str(a),str(b)],area=enc(A),physical=True,unperturbed_upper=enc(-A)));continue
        deficit=(target-A).upper()
        assert gamma*deficit<q('.3')
        necessary=2*deficit if free else deficit/(1-gamma*deficit)
        H3=[unit(3*t) for t in T];H5=[unit(5*t) for t in T]
        M=[];seen=set()
        for i,j,fix,mask in patterns([float(x) for x in mids],[float(x.mid()) for x in L]):
            w=[arb(0)]*8
            for k,bit in zip(fix,mask):w[k]=L[k]*int(bit)
            rhs=tuple(-x for x in vsum(N,w));det=cross(N[i],N[j])
            w[i]=cross(rhs,N[j])/det;w[j]=cross(N[i],rhs)/det
            w=[(1-blend)*x+blend*l/2 for x,l in zip(w,L)]
            if not all(x>0 and l-x>0 for x,l in zip(w,L)):continue
            A3=vsum(H3,w);A5=vsum(H5,w)
            key=tuple(round(float(x.mid()),6) for x in (*A3,*A5))
            if key in seen:continue
            seen.add(key); M.append((*A3,*A5))
        assert len(M)>12
        Vf=np.array([[float(x.mid()) for x in row] for row in M])
        # For large negative skew the pure cubic has larger clearance; in all
        # nonnegative skew cells a single fixed two-harmonic body is used.
        model='cubic' if b<=F('-.85') else 'mixed'
        w3=s/16 if model=='cubic' else s*q('7/40')/16
        w5=arb(0) if model=='cubic' else s*q('33/40')/48
        best=arb(1);leaf=0;minwidth=F(360)
        for hand in ([1] if not w5 else [1,-1]):
            # With reflection, only the fixed relative phase changes sign.
            todo=[(F(j),F(j+10)) for j in range(0,360,10)]
            while todo:
                x,y=todo.pop();mid=(x+y)/2
                angle=float(mid)*pi/180;hh=hand*144*pi/180
                proj=np.array([float(w3)*cos(3*angle),float(w3)*sin(3*angle),
                               float(w5)*cos(5*angle+hh),float(w5)*sin(5*angle+hh)])
                z=M[int(np.argmax(Vf@proj))]
                u=deg(mid);delta=deg((y-x)/2)
                t3=3*u;t5=5*u+hand*phase
                f=w3*(z[0]*t3.cos()+z[1]*t3.sin())+w5*(z[2]*t5.cos()+z[3]*t5.sin())
                df=3*w3*(-z[0]*t3.sin()+z[1]*t3.cos())+5*w5*(-z[2]*t5.sin()+z[3]*t5.cos())
                d2=9*w3*(abs(z[0])+abs(z[1]))+25*w5*(abs(z[2])+abs(z[3]))
                score=(f-abs(df)*delta-d2*delta**2/2).lower()
                if score>necessary:
                    best=min(best,score);leaf+=1;minwidth=min(minwidth,y-x);continue
                assert y-x>F('.008'),(panel,a,b,model,hand,x,y,score,necessary,len(M))
                todo.extend(((x,mid),(mid,y)))
        if not free:
            best=min(best,(q('.9')/gamma).lower())
            assert gamma*best<1
        wall=A+(best/2 if free else best/(1+gamma*best))
        assert wall>target,(a,wall,target)
        results.append(dict(panel=panel,d=[str(a),str(b)],area=enc(wall),bare_area=enc(A),unperturbed_upper=enc(-A),
                    body=model,angle_leaves=leaf,minwidth=str(minwidth),dual_count=len(M),dual_floor=enc(best)))
    return dict(complete=True,d=[str(lo),str(hi)],normal_errors=str(er),arbitrary_convex_extension=free,
                projection_bounds='0 <= weights <= facet lengths; exact vector equilibrium',harmonic_factor='.999',mix=['7/40','33/40'],
                relative_phase_deg=144,area_floor=str(min(F(row['area']) for row in results)),
                target=enc(target),threshold_B=enc(B),panels=results,python=sys.version,flint=flint_version,numpy=np.__version__,
                script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),precision=ctx.prec)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for key in ('low','high','error','target'):ap.add_argument(key,type=F)
    ap.add_argument('panels',type=int);ap.add_argument('output',type=Path);ap.add_argument('--free',action='store_true')
    a=ap.parse_args();r=run(a.low,a.high,a.error,a.panels,a.free,a.target)
    a.output.write_text(json.dumps(r,indent=1)+'\n');print(float(F(r['area_floor'])),sum(z.get('angle_leaves',0) for z in r['panels']))
