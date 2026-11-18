"""Exact global certificate for C4 augmented with the complete Reuleaux pentagon."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
import flint
from flint import arb, ctx
ctx.prec=192
pi=arb.pi()
half=arb(1)/2
root=arb(3).sqrt()
ell=half/root
inner=pi/8+root/4

def q(x):
    x=F(x)
    return arb(x.numerator)/x.denominator
def rad(x):
    return q(x)*pi/180
def unit(x):
    return x.cos(),x.sin()
def vec(k,v):
    return k*v[0],k*v[1]
def add(a,b):
    return a[0]+b[0],a[1]+b[1]
def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]
def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]
def rotate(s,v):
    z=unit(s)
    return dot((z[0],-z[1]),v),dot((z[1],z[0]),v)
def moment(l,h):
    return h.sin()-l.sin(),l.cos()-h.cos()
def bounds(x):
    def asrat(y):
        m,k=y.man_exp()
        return str(F(int(m))*F(2)**int(k))
    return [asrat(x.lower()),asrat(x.upper())]
r=half/rad(18).cos()
us=[unit(rad(60*i)) for i in range(3)]

def make_poly(c,C,R,wall):
    angle=rad(c)
    verts=[vec(r,unit(angle+rad(72*j))) for j in range(5)]
    V=(arb(0),arb(0)); H=(arb(0),arb(0))
    # Fixed host blocks give C0 + C1 cos s + C2 sin s.
    C0=arb(0);C1=arb(0);C2=arb(0)
    for l,h,j,a in C:
        v=vec(half,moment(l,h)); mass=(h-l)/2
        V=add(V,v); C0-=half*mass
        C1+=dot(verts[j],v)+a*mass
        C2+=cross(verts[j],v)
    for l,h,j,a in R:
        v=moment(l,h);H=add(H,v);C0+=dot(verts[j],v)+half*(h-l)
    i,sgn,atten=wall; k,m=[z for z in range(3) if z!=i]
    atten=q(atten);sgn=q(sgn)
    def ev(s):
        co,si=unit(s)
        G=vec(atten,add(V,rotate(s,H)))
        tot=atten*(C0+C1*co+C2*si)
        w=[None,None,None];w[i]=ell*sgn
        tar=vec(-1,add(G,vec(w[i],us[i])))
        det=cross(us[k],us[m])
        w[k]=cross(tar,us[m])/det;w[m]=cross(us[k],tar)/det
        b=[r*(angle+s).cos()-half,half+r*(angle+s+rad(156)).cos(),
           half+r*(angle+s+rad(168)).cos()]
        return inner+(tot+sum(x*y for x,y in zip(w,b)))/2,w[k],w[m]
    values=[ev(rad(t))[0] for t in (0,90,180,270,45)]
    A=sum(values[:4])/4
    B=(values[0]-values[2])/2;C=(values[1]-values[3])/2
    D=(values[0]+values[2]-values[1]-values[3])/4
    E=values[4]-A-(B+C)/arb(2).sqrt()
    def second(s):
        return -B*s.cos()-C*s.sin()-4*(D*(2*s).cos()+E*(2*s).sin())
    return ev,second

def replay(src):
    x=r*rad(6).cos()-half
    p=(-x,x/root)
    qr=[add(vec(r,unit(rad(6+72*j))),p) for j in range(5)]
    b=(half/dot(qr[4],qr[4]).sqrt()).acos()
    U=rad(120)-2*b;L=rad(300)-U
    assert rad(104)<U<rad(105)
    assert cross(vec(-1,qr[4]),unit(rad(120)))>0
    assert dot(qr[1],unit(rad(96)))>half
    assert cross(qr[4],unit(rad(300)))>0 and dot(qr[4],unit(rad(312)))>0
    assert (r+2*x/root)*rad(18).cos()>half
    i1=dot(qr[1],moment(rad(60),rad(96)))-half*rad(36)
    i2=half*(U-rad(96))+dot(qr[4],moment(rad(96),U))
    target=inner+(i1+3*i2)/2
    c6=[(rad(60),rad(96),1,0),(rad(96),U,4,1),(L,rad(204),0,1),(rad(204),rad(240),3,0)]
    r6=[c6[1],c6[2]]
    specs=[(F(6),F('6.05'),F(6),c6,r6,(0,1,1))]
    def inside(l,h,centers,hop):
        return any(v-hop<=l<=h<=v+hop for v in centers)
    for item in json.loads(src.read_text()):
        lo,hi=map(F,item['interval']);c=(lo+hi)/2
        assert lo==specs[-1][1] and lo<hi
        blocks=[]
        for kind in ('C','R'):
            group=[];prev=F(-1)
            for x1,x2,j,a in item[kind]:
                l=F(x1,10**7);h=F(x2,10**7)
                assert prev<=l<h and j in range(5) and a in (0,1)
                prev=h
                if a:
                    assert inside(l,h,[c+72*j+180+360*k for k in range(-3,3)],18)
                if kind=='C':
                    assert inside(l,h,(90,210,330),30)
                else:
                    assert a==1 and inside(l+lo-c,h+hi-c,(90,210,330),30)
                group.append((rad(l),rad(h),j,a))
            blocks.append(group)
        i,sgn,atten=item['wall']
        assert i in range(3) and -1<=F(sgn)<=1 and 0<=F(atten)<=1
        specs.append((lo,hi,c,*blocks,(i,sgn,atten)))
    assert specs[-1][1]==18
    checks=[]
    for idx,(lo,hi,c,C,R,wall) in enumerate(specs):
        ff,sec=make_poly(c,C,R,wall)
        if idx==0:
            assert ff(rad(0))[0].overlaps(target)
        l,h=rad(lo-c),rad(hi-c)
        box=l.union(h)
        _,v,w=ff(box)
        gap=min((ell-abs(v)).lower(),(ell-abs(w)).lower())
        deriv=sec(box).upper()
        assert gap>0 and deriv<0,(lo,hi,gap,deriv)
        left=ff(l)[0]-target;right=ff(h)[0]-target
        assert right>0 and (idx==0 or left>0)
        checks.append(dict(interval=[str(lo),str(hi)],left_endpoint=bounds(left) if idx else 'exact contact',
                           right_endpoint=bounds(right),free_wall_slack=bounds(gap),
                           second_derivative_upper=bounds(deriv)))
    # Separating direction in the previously certified optimal five-vertex hull.
    p18=2*x/root
    Q64=(r+p18)*rad(64).sin()
    R64=1+r*rad(170).cos()+p18*rad(64).sin()
    assert Q64<half and R64-half>q('.00597166')
    qa=(half/(r+p18)).acos()
    vertex3=add(vec(r,unit(rad(234))),(arb(0),p18))
    qb=(half/dot(vertex3,vertex3).sqrt()).acos()
    pentagon_only=inner+(qa.tan()-qa+2*(qb.tan()-qb))/4
    difference=target-pentagon_only
    print('full minimum',target,'increase from five-vertex minimum',difference)
    print('closed angle panels',len(specs))
    return dict(full_minimum=bounds(target),translation=[bounds(z) for z in p],central_breakpoint_degrees=bounds(U*180/pi),
         increase_from_five_vertices=bounds(difference),exposed_at_64_deg=bounds(R64-half),five_vertex_support_at_64_deg=bounds(Q64),panels=checks,
         precision=ctx.prec,python=sys.version,flint=flint.__version__,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
         panel_sha256=sha256(src.read_bytes()).hexdigest())
if __name__=='__main__':
    result=replay(Path(__file__).with_name('panels.json'))
    Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
