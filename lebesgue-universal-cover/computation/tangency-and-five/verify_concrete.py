"""Exact fixed-host miss of the regular pentagon and support forced by a fifth harmonic."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import itertools,json,sys
import numpy as np
from flint import arb,ctx,__version__ as flint_version
ctx.prec=192
def Q(v):
    v=F(v);return arb(v.numerator)/v.denominator
def deg(v):return Q(v)*arb.pi()/180
def u(t):return t.cos(),t.sin()
def uv(v):return u(deg(v))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def rot(a):return -a[1],a[0]
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(k,a):return k*a[0],k*a[1]
def norm(a):return dot(a,a).sqrt()
def balance(a,b,c):
    v=[cross(b,c),cross(c,a),cross(a,b)]
    assert min(v)>0
    return [x/sum(v) for x in v]
def bound(q):
    ret=[]
    for x in (q.lower(),q.upper()):
        m,k=x.man_exp();ret.append(str(F(int(m))*F(2)**int(k)))
    return ret
def run():
    r=1/(2*deg(18).cos());root=arb(3).sqrt();half=Q('.5');e6=uv(6);v330=uv(330)
    b=r*(root*e6[0]/2+deg(72).cos());d=r**2*e6[1]**2
    s=d/(b+(b**2-d/4).sqrt())
    delta=r*e6[0]+root*s/2-half
    assert 0<s<delta
    t=mul(s,v330);ys=[add(mul(r,uv(6+72*j)),t) for j in range(5)]
    assert 0<ys[0][1]<1/(2*root)
    assert 0<r-s<1/root
    # The contacts of vertices 1 and 3 are interior to two surviving disk arcs.
    v1=mul(1/(half+delta),ys[1]);v3=mul(1/(half+delta),ys[3])
    for p,l,h in ((v1,60,120),(v3,180,240)):
        assert cross(uv(l),p)>0 and cross(p,uv(h))>0
    v4=uv(300)
    weights=balance(v1,v3,v4);dirs=(v1,v3,v4);ids=(1,3,4)
    sector=deg(6).union(deg(F('6.25')))
    derivative=sum(w*r*dot(n,rot(u(sector+deg(72*k)))) for w,n,k in zip(weights,dirs,ids))
    assert derivative>0
    # Beyond the exact endpoint use only fixed rational-degree support normals.
    degrees=[i for i in range(0,360,6) if i%120==0 or i%120>=60]
    dirs=[uv(z) for z in degrees]; sf=np.deg2rad(np.array(degrees)); R=float(r)
    choices=[];weight_values=[]
    for idx in itertools.combinations(range(len(degrees)),3):
        v=[dirs[i] for i in idx]
        if min(cross(v[0],v[1]),cross(v[1],v[2]),cross(v[2],v[0]))<Q('.001'):continue
        w=balance(*v)
        weight_values.append([float(x) for x in w]);choices.append(idx)
    idxf=np.array(choices); wf=np.array(weight_values)
    cells=[]; margins=[]
    for j in range(47):
        lo=(F(25)+j)/4;hi=lo+F(1,4);mid=(lo+hi)/2
        phis=2*np.pi*np.arange(5)/5+np.deg2rad(float(mid))
        vertex=np.argmax(np.cos(sf[:,None]-phis[None,:]),axis=1)
        hs=R*np.cos(sf-phis[vertex])-.5
        k=int(np.argmax(np.sum(wf*hs[idxf],axis=1)))
        ix=choices[k];angles=[degrees[i] for i in ix];ks=[int(vertex[i]) for i in ix];w=balance(*[dirs[i] for i in ix])
        theta=deg(lo).union(deg(hi))
        assert all(dot(dirs[i],u(theta+deg(72*z)))>0 for i,z in zip(ix,ks))
        f=lambda x:sum(co*r*dot(dirs[i],uv(x+72*z)) for co,i,z in zip(w,ix,ks))-half
        margin=min((f(lo)-delta).lower(),(f(hi)-delta).lower())
        assert margin>0,(lo,hi,angles,ks,margin)
        cells.append(dict(degrees=[str(lo),str(hi)],normals=angles,vertices=ks,margin=bound(margin)[0]));margins.append(margin)
    assert F(cells[-1]['degrees'][1])==18
    c=Q('999/48000')
    h1=half+root*c*deg(F('4.62')).cos()
    h10=half+2*c*deg(F('34.62')).cos()*deg(50).cos()
    h15=half+2*c*deg(F('49.62')).cos()*deg(60).cos()
    return dict(complete=True,intrinsic_radius=bound(r),translation_radius=bound(s),
                exact_host_miss=bound(delta),endpoint_derivative_per_radian=bound(derivative),
                remaining_angle_margin=bound(min(margins)),angle_panels=cells,
                gibbs_sigma_degrees=['1.50','1.54'],fifth_harmonic_width_one_amplitude='999/48000',
                third_cap_lower=bound(h1),whole_10_degree_arc_lower=bound(h10),whole_15_degree_arc_lower=bound(h15),
                python=sys.version,flint=flint_version,numpy=np.__version__,precision=ctx.prec,script_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    j=run();Path(sys.argv[1]).write_text(json.dumps(j,indent=1)+'\n')
    print({k:float(F(j[k][0])) for k in ('exact_host_miss','translation_radius','endpoint_derivative_per_radian','remaining_angle_margin','third_cap_lower','whole_10_degree_arc_lower','whole_15_degree_arc_lower')})
