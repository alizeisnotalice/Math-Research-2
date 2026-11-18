"""Recompute certified hulls and translation obstructions for all angles and both hands."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import gzip
import json
import sys
import flint
from flint import arb,ctx
sys.path.insert(0,str(Path(__file__).resolve().parent))
from template import build

def q(x):
    x=F(x)
    return arb(x.numerator)/x.denominator

def det(x,y):
    return x[0]*y[1]-x[1]*y[0]

def balance(ns):
    if len(ns)==3:
        w=[det(ns[(i+1)%3],ns[(i+2)%3]) for i in range(3)]
        if w[0]<0:w=[-a for a in w]
    else:
        assert len(ns)==2
        k=0 if ns[0][0] else 1
        w=[abs(ns[1][k]),abs(ns[0][k])]
    assert min(w)>0 and all(sum(n[k]*z for n,z in zip(ns,w))==0 for k in range(2))
    return w

def read(path):
    raw=path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)

def verify(c,source,data):
    ctx.prec=160
    raw=source.read_bytes()
    assert sha256(raw).hexdigest()=='45cb554f79d2c30f23cdc526f7cc1a333b6f98c3dbaa60234862d449bb66de05'
    cloud=json.loads(raw)
    P=[tuple(F(z,cloud['denominator']) for z in row) for row in cloud['points']]
    diameter=max((x-X)**2+(y-Y)**2 for x,y in P for X,Y in P)
    assert len(P)==11 and diameter<=1 and max(x*x+y*y for x,y in P)<F('.58')**2
    points=[(q(x),q(y)) for x,y in P]
    paired='double_cases' in data
    h,D=build(c,data.get('seven',False),data['samples'],data.get('five_only',paired),box=data.get('box'))
    n=len(h)
    pre=[0]
    for i in range(2*n):pre.append(pre[-1]+det(h[i%n],h[(i+1)%n]))
    assert pre[n]>0
    A=F(data['target'])
    E=2*A*D**2-pre[n]
    if data['cells']==0:
        assert E<0
        return dict(target=str(A),direct_polygon_area=str(F(pre[n],2*D*D)),cases=0)
    assert E>=0
    R=q('.706');radius=q('.5')
    single_point=(arb.pi()-(radius/R).acos())*radius**2+radius*(R**2-radius**2).sqrt()
    if paired:assert single_point>q(A)
    @lru_cache(None)
    def chord(i,j):
        assert 0<=i<n and 0<=j<n and i!=j
        x,y=h[i];X,Y=h[j];u,v=Y-y,x-X
        length=arb(u*u+v*v).sqrt()
        return (u,v),length,pre[j if j>i else j+n]-pre[i],q(F(u*x+v*y,D))/length
    @lru_cache(None)
    def disjoint(a,b):
        (u,v),la,_,ba=chord(*a);(U,V),lb,_,bb=chord(*b)
        return ba+bb>R*((u/la+U/lb)**2+(v/la+V/lb)**2).sqrt()
    K=data['period'];ends={-1:0,1:0};gap=None;doubles=0
    for sign,lo,hi,systems in data['cells']:
        assert sign in ends and lo==ends[sign] and lo<hi<=K
        ends[sign]=hi
        if not paired:systems=[[z] for z in systems]
        ds=any(len(z)==2 for z in systems);doubles+=ds
        angle=arb.pi()*q(F(lo+hi,K));sn,cs=angle.sin_cos()
        loss=q(F('.58')*(hi-lo)/K)*arb.pi()
        normals=[];rhs=[];scales=[]
        for group in systems:
            assert 1<=len(group)<=2
            if len(group)==2:assert disjoint(tuple(group[0]),tuple(group[1]))
            norm=[0,0];scale=arb(0);upper=E;support=arb(0)
            for i,j in group:
                (u,v),length,arc,_=chord(i,j)
                xx=cs*u+sn*v;yy=sign*(-sn*u+cs*v)
                support+=max((xx*px+yy*py).lower() for px,py in points)
                norm[0]+=u;norm[1]+=v;scale+=length;upper+=arc
            normals.append(norm);scales.append(scale);rhs.append(support-q(upper/D)-loss*scale)
        weights=balance(normals)
        val=sum(w*r for w,r in zip(weights,rhs))
        assert val>0,(sign,lo,hi,val)
        margin=(val/sum(w*scale for w,scale in zip(weights,scales))).lower()
        gap=margin if gap is None else min(gap,margin)
    assert all(z==K for z in ends.values()) and doubles==data.get('double_cases',0)
    budgets={}
    for target in map(F,('.83750','.83755')):
        if A<=target:continue
        a,b=0,10**8
        while b-a>1:
            m=(a+b)//2
            if A>target*(1+F(2*m,10**8))**2:a=m
            else:b=m
        budgets[str(target)]=str(F(a,10**8))
    return dict(target=str(A),minimum_normalized_margin=str(gap),cases=len(data['cells']),
                disjoint_cap_cases=doubles,independent_motion_budgets=budgets,
                diameter_squared=str(diameter),inner_polygon_area=str(F(pre[n],2*D*D)),hull_vertices=n)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('computation',type=Path);p.add_argument('folder',type=Path)
    p.add_argument('--output',type=Path)
    a=p.parse_args();assert __debug__
    root=a.folder;c=a.computation
    result={}
    for f in sorted(root.glob('*.json.gz')):
        data=read(f)
        if f.name.startswith('star-'):
            sys.path.insert(0,str(c.resolve()))
            import certify_five_body_contact_stationary as st
            ctx.prec=192
            _,actual,_,_=st.model(*(st.interval(*v) for v in st.BOX))
            assert all(q(F(lo,2**44))<=x<=q(F(hi,2**44)) for (lo,hi),x in zip(data['box'],actual))
        result[f.name]=verify(c,root/'new-points.json',data)
        print(f.name,result[f.name],flush=True)
    assert result
    if a.output:a.output.write_text(json.dumps(dict(python_flint=flint.__version__,cases=result),indent=2)+'\n')
