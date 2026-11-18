"""Produce and certify two 60-point rational subsets of the exact inputs."""
from pathlib import Path
import sys,json,hashlib,platform
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from cross_common import *

OUT=Path(__file__).resolve().parent
S=10**12
shrink=F(999999999,1000000000)
vertices,_=body(ns,F('.7982907'),F('.032688'))
edges=[sub(vertices[i],vertices[j]) for i,j in ns['orders']]

def pent_contact(u):
    cr=[ns['cross'](e,u) for e in edges]
    good=[j for j in range(10) if cr[j]>0 and cr[(j+1)%10]<0]
    assert len(good)==1, ('pent sector',good)
    k,r=ns['sector'][good[0]]
    return add(vertices[k],mul(u,arb(r)))

def five_contact(t):
    # Use exact antipodal correspondence on the second semicircle.
    i=0 if t<pi else 1
    a=t-i*pi
    eligible=[j for j in range(6) if a>xs[j] and a<xs[j+1]]
    assert len(eligible)==1, ('five sector',t,eligible)
    r,c=pieces[eligible[0]];u=unit(a)
    x=add(c,mul(u,arb(r)))
    return sub(x,u) if i else x

def rounding(x):
    y=q(shrink)*x
    k=round(float(y)*S)
    assert abs(y*S-k)<1
    return k

# Prove .42 B is contained, by 128 equally spaced support directions and
# the .58 radius established by cross_common. Every direction is <= pi/128
# from a checked normal, and support variation is <= .58*pi/128.
for name,fn in [('r151',hp),('r216',hf)]:
    vals=[fn(unit(2*pi*i/128+q(F(371,1000)))) for i in range(128)]
    assert all(h > A('.42')+A('.58')*pi/128 for h in vals),name

assert q(F(1,10**9))*A('.42') > A('1.5')/S
out={}
for name in ('r151','r216'):
    points=[];errors=[]
    for i in range(60):
        t=q(F(371,1000))+2*pi*i/60
        u=unit(t)
        # Reduce only the angle passed to five_contact, not the normal.
        tt=t if t<2*pi else t-2*pi
        v=pent_contact(u) if name=='r151' else five_contact(tt)
        fn=hp if name=='r151' else hf
        # An explicit independent support check for the selected contact.
        assert (dot(v,u)-fn(u)).contains(0),(name,i,dot(v,u),fn(u))
        z=[rounding(v[0]),rounding(v[1])]
        assert all(abs(q(shrink)*v[k]-q(F(z[k],S)))<q(F(1,S)) for k in (0,1))
        points.append(z)
    # Repeated support normals may touch the same corner. Retain one copy.
    points=list(dict.fromkeys(map(tuple,points)))
    perimeter=sum(arb((a[0]-b[0])**2+(a[1]-b[1])**2).sqrt()/S for a,b in zip(points,points[1:]+points[:1]))
    assert perimeter < pi
    d2=max((x[0]-y[0])**2+(x[1]-y[1])**2 for x in points for y in points)
    assert d2<=S*S,(name,d2)
    out[name]={'denominator':S,'points':points,'diameter_squared_numerator':d2,'cyclic_tour_length_upper':str(perimeter.upper())}
    print(name,'points',len(points),'tour',perimeter.upper(),'diameter2',d2,'/',S*S,flush=True)
obj={'offset_radians':'371/1000','count':60,'shrink':str(shrink),'inscribed_radius':'21/50','coordinate_error_bound':'1/1000000000000','sets':out,'flint':flint.__version__,'precision_bits':ctx.prec}
target=OUT/'points60.json'
if target.exists():
    saved=json.loads(target.read_text())
    assert saved['sets']==json.loads(json.dumps(obj['sets'])) and saved['offset_radians']==obj['offset_radians'] and saved['shrink']==obj['shrink']
else:
    target.write_text(json.dumps(obj,indent=2)+'\n')
