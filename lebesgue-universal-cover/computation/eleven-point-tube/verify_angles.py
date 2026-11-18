"""Prove uniform angle/translation area and continuous host-pose perturbation bounds."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
import flint
from flint import arb,ctx

sys.path.insert(0,str(Path(__file__).resolve().parent))
from template import build


def rational(x):
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


def verify(c,p,data):
    ctx.prec=160
    assert sha256(p.read_bytes()).hexdigest()=='05fb76a261b46c7b939885becb32b31b3e5faca2b017dde5ff70123ff22f6a0d'
    raw_points=json.loads(p.read_bytes())['points']
    xyz=[tuple(map(F,z)) for z in raw_points]
    assert len(xyz)==11 and max(x*x+y*y for x,y in xyz)<F('.58')**2
    diameter=max((x-X)**2+(y-Y)**2 for x,y in xyz for X,Y in xyz)
    assert diameter<=1
    points=[tuple(map(rational,z)) for z in xyz]
    h,D=build(c,data['seven'],data['samples'],data.get('five_only',False),data.get('omit',[]),data.get('fraction'))
    n=len(h);twice=sum(det(h[i],h[(i+1)%n]) for i in range(n));assert twice>0
    cross=[0]
    for i in range(2*n):cross.append(cross[-1]+det(h[i%n],h[(i+1)%n]))
    K=data['period'];A=F(data['target']);E=2*A*D**2-twice
    assert E>0
    ends={-1:0,1:0};gap=None;max_levels=0
    for sign,lo,hi,arcs in data['cells']:
        assert sign in ends and lo==ends[sign] and lo<hi<=K
        ends[sign]=hi
        angle=arb.pi()*rational(F(lo+hi,K));sn,cs=angle.sin_cos()
        loss=rational(F('.58')*(hi-lo)/K)*arb.pi()
        normals=[];rhs=[]
        for i,j in arcs:
            assert 0<=i<n and 0<=j<n and i!=j
            x,y=h[i];X,Y=h[j]
            normals.append((Y-y,x-X))
            rhs.append((E+cross[j if j>i else j+n]-cross[i])/D)
        weights=balance(normals);value=normsum=arb(0)
        for (u,v),upper,w in zip(normals,rhs,weights):
            x=cs*u+sn*v;y=sign*(-sn*u+cs*v)
            support=max((x*px+y*py).lower() for px,py in points)
            norm=arb(u*u+v*v).sqrt()
            value+=w*(support-rational(upper)-loss*norm)
            normsum+=w*norm
        assert value>0,(sign,lo,hi,value)
        margin=(value/normsum).lower()
        gap=margin if gap is None else min(gap,margin)
    assert all(z==K for z in ends.values())
    # Each B_n has constant width 1 and circumradius 1/(2*cos(pi/(2*n))).
    # Thus |h'| <= (1/2)*tan(pi/(2*n)), giving legal translation + angle budgets.
    L={3:'.289',5:'.163',7:'.115',9:'.089'}
    assert all((arb.pi()/(2*k)).tan()/2<rational(v) for k,v in L.items())
    sample=json.loads((c/'finite-point-leaf-input.json').read_bytes())
    pp=[(F(x,sample['denominator']),F(y,sample['denominator'])) for x,y in sample['points']]
    assert all((x+F('.01007'))**2+y*y<F('.529')**2 for x,y in pp)
    sample=json.loads((c/'multi-host-points.json').read_bytes())
    qq=[(F(x,sample['point_denominator']),F(y,sample['point_denominator'])) for x,y in sample['points']]
    cx,cy=[sum(z[i] for z in qq)/len(qq) for i in (0,1)]
    qr=max((x-cx)**2+(y-cy)**2 for x,y in qq)
    assert qr<F('.60')**2
    tr={}
    for target in (F('.83750'),F('.83755')):
        if A<=target:continue
        # Compute a certified 1e-8-rounded safe radius by integer comparisons alone.
        e=0;step=1
        while A>target*(1+F(2*step,10**8))**2:step*=2
        while step:
            if A>target*(1+F(2*(e+step),10**8))**2:e+=step
            step//=2
        eps=F(e,10**8)
        assert A/(1+2*eps)**2>target
        tr[str(target)]=dict(motion_budget=str(eps),lower=str(A/(1+2*eps)**2))
    return dict(model='fixed disk; other participating bodies may move independently',
                seven=data['seven'],fraction=data.get('fraction'),omit=data.get('omit',[]),five_only=data.get('five_only',False),samples=data['samples'],
                host_polygon_vertices=n,host_inner_area=str(F(twice,2*D**2)),
                all_placed_point_enlarged_area_strict_lower=str(A),diameter_squared=str(diameter),
                cases=len(data['cells']),minimum_normalized_strict_margin=str(gap),motion_budgets=tr,
                regular_angular_constants=L,point_radius_about_minus_01007='.529',
                Q_centered_radius_squared=str(qr),points_sha256=sha256(p.read_bytes()).hexdigest(),
                arb_precision=ctx.prec,python_flint=flint.__version__)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('computation','points','certificate','output'):parser.add_argument(name,type=Path)
    a=parser.parse_args();assert __debug__
    result=verify(a.computation,a.points,json.loads(a.certificate.read_bytes()))
    result['certificate_sha256']=sha256(a.certificate.read_bytes()).hexdigest()
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['all_placed_point_enlarged_area_strict_lower'],result['motion_budgets'])
