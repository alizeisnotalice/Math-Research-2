"""Explicit two-parameter contact system near the four-body feasible minimum."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb,arb_series,ctx
from certify_three_body_symmetric_slice import add,det,dot,interval,pi,q,r,scale,sub

ctx.cap=3


def left(x):
    return -x[1],x[0]


def normal(t):
    return t.cos(),t.sin()


def lin2(a,b,u,v):
    d=det(a,b)
    return (u*b[1]-v*a[1])/d,(v*a[0]-u*b[0])/d


def unit_root(p,h,sign):
    modulus=dot(p,p)
    return add(scale(h/modulus,p),scale(sign*(modulus-h*h).sqrt()/modulus,left(p)))


def other(p,h,n):
    return sub(scale(2*h/dot(p,p),p),n)


def regular(n,k):
    return scale(1/(2*(pi/(2*n)).cos()),normal(2*pi*k/n))


def motion(n,a,b,c,specs):
    index0,level0=specs[0]
    index1,level1=specs[1]
    index2,level2=specs[2]
    v0,v1,v2=(regular(n,k) for k in (index0,index1,index2))
    t0=lin2(a,b,level0,level1)
    tc=lin2(a,b,-dot(a,v0),-dot(b,v1))
    ts=lin2(a,b,-dot(a,left(v0)),-dot(b,left(v1)))
    linear=(dot(c,add(tc,v2)),dot(c,add(ts,left(v2))))
    rhs=level2-dot(c,t0)
    cs,sn=unit_root(linear,rhs,-1)
    if isinstance(cs,arb):
        assert cs>0
    t=add(t0,add(scale(cs,tc),scale(sn,ts)))
    theta=(sn/cs).atan()
    return theta,t, [add(t,add(scale(cs,v),scale(sn,left(v)))) for v in [regular(n,k) for k in range(n)]]


def geometry(alpha,beta):
    A,B=normal(alpha),normal(beta)
    t3=lin2(A,B,r-dot(A,regular(3,0)),r-dot(B,regular(3,1)))
    p3=[add(t3,regular(3,k)) for k in range(3)]
    D=other(p3[0],r,A)
    theta7,t7,p7=motion(7,A,B,D,[(0,r),(1,r),(2,-r)])
    C=other(p7[2],-r,D)
    theta5,t5,p5=motion(5,A,B,C,[(3,-r),(1,r),(2,-r)])
    placements=(*t3,theta5,*t5,theta7,*t7)
    normals={'A':A,'B':B,'C':C,'D':D,'E':other(p7[0],r,A),'F':other(p5[3],-r,A),
             'X':unit_root(sub(p5[2],p3[0]),1,-1),'Y':other(p5[2],r,scale(-1,C)),
             'Z':unit_root(p7[3],r,-1)}
    return placements,normals,{3:p3,5:p5,7:p7}


def ratio_angle(n):
    return (n[1]/n[0]).atan()


def contacts(alpha,beta,place,ns):
    return {'A':alpha,'B':beta,'C':2*pi+ratio_angle(ns['C']),'D':2*pi+ratio_angle(ns['D']),
        'E':ratio_angle(ns['E']),'F':ratio_angle(ns['F']),
        **{name: pi+ratio_angle(ns[name]) for name in ['X','Y','Z']}}


def chain(alpha,beta):
    placements,ns,p=geometry(alpha,beta)
    ev=contacts(alpha,beta,placements,ns)
    phi5,phi7=placements[2],placements[5]
    origin=(0*alpha,0*alpha)
    pieces=[
        ('A',ev['A'],3,0,p[3][0]),
        ('E',ev['E'],7,0,p[7][0]),
        ('F',ev['F'],0,r,origin),
        ((5,1,-1),3*pi/10+phi5,5,1,p[5][3]),
        ('B',ev['B'],5,0,p[5][1]),
        ((3,1,1),5*pi/6,3,0,p[3][1]),
        ('X',ev['X'],3,1,p[3][0]),
        ('Y',ev['Y'],5,0,p[5][2]),
        ('Z',ev['Z'],0,r,origin),
        ((7,3,1),pi+phi7-pi/14,7,0,p[7][3]),
        ('Aopp',pi+alpha,7,1,p[7][0]),
        ((3,2,-1),7*pi/6,3,1,p[3][0]),
        ((3,2,1),3*pi/2,3,0,p[3][2]),
        ('Bopp',pi+beta,3,1,p[3][1]),
        ((5,4,-1),3*pi/2+phi5,5,1,p[5][1]),
        ((5,4,1),17*pi/10+phi5,5,0,p[5][4]),
        ('C',ev['C'],5,1,p[5][2]),
        ('D',ev['D'],0,r,origin),
        ('end',2*pi,3,0,p[3][0]),
    ]
    return pieces,placements,ns,p


def area(alpha,beta):
    pieces,_,_,_=chain(alpha,beta)
    start=0*alpha; value=0*alpha
    for (_,stop,_,rr,c),(_,_,_,rr_next,cn) in zip(pieces,pieces[1:]+pieces[:1]):
        n0,n1=normal(start),normal(stop)
        value+=rr*det(c,sub(n1,n0))+rr**2*(stop-start)
        value+=det(add(c,scale(rr,n1)),add(cn,scale(rr_next,n1)))
        start=stop
    return value/2


def jets(a,b):
    x=area(arb_series([a,1]),arb_series([b]))
    y=area(arb_series([a]),arb_series([b,1]))
    z=area(arb_series([a,1]),arb_series([b,1]))
    hxy=z[2]-x[2]-y[2]
    return (x[1],y[1]),((2*x[2],hxy),(hxy,2*y[2]))



BOX=(('.51250444532618500354390','.51250444532618500354394'),
     ('1.57747301538059498239155','1.57747301538059498239160'))
PIECE_DATA=((3,0,0),(7,0,0),(0,0,r),(5,3,1),(5,1,0),(3,1,0),(3,0,1),
            (5,2,0),(0,0,r),(7,3,0),(7,0,1),(3,0,1),(3,2,0),(3,1,1),
            (5,1,1),(5,4,0),(5,2,1),(0,0,r),(3,0,0))
TIES={name:{0,3,5,7} for name in ('A','Aopp','B','Bopp')}
TIES.update({'C':{0,5,7},'D':{0,3,7},'E':{0,7},'F':{0,5},'X':{3,5},'Y':{0,5},'Z':{0,7}})


def primitive_key(n,angle,rotation):
    turn=2*pi/n
    j=round(float((angle-rotation)/turn))
    gap=angle-rotation-j*turn
    if abs(gap)<pi/(2*n):
        return n,j%n,0
    j=round(float((angle-rotation-pi)/turn))
    assert abs(angle-rotation-pi-j*turn)<pi/(2*n)
    return n,j%n,1


def verify_topology(a,b):
    pieces,t,ns,p=chain(a,b)
    assert all(ns[k][0]>0 for k in ('C','D','E','F'))
    assert all(ns[k][0]<0 for k in ('X','Y','Z'))
    cuts={'start':arb(0),'end':2*pi}
    rotations={3:arb(0),5:t[2],7:t[5]}
    for n,f in rotations.items():
        for i in range(n):
            for end in (-1,1):
                angle=f+2*pi*i/n+end*pi/(2*n)
                if angle<0: angle+=2*pi
                if angle>2*pi: angle-=2*pi
                assert 0<angle<2*pi
                cuts[n,i,end]=angle
    for k,angle,*_ in pieces:
        if k not in cuts:
            cuts[k]=angle
    keys=sorted(cuts,key=lambda k:float(cuts[k].mid()))
    assert [*keys[:1],*keys[-1:]]==['start','end']
    assert all(cuts[k]<cuts[j] for k,j in zip(keys,keys[1:]))
    active=0
    endpoints=0
    critical=0
    for k,j in zip(keys,keys[1:]):
        l,h=cuts[k],cuts[j]
        assert h-l<pi
        data=PIECE_DATA[active]
        midpoint=(l+h)/2
        if data[0] in rotations:
            assert primitive_key(data[0],midpoint,rotations[data[0]])==data
        owner,vertex,rad=data
        center=p[owner][vertex] if owner else (arb(0),arb(0))
        nl,nh=normal(l),normal(h)
        for n in (0,3,5,7):
            if owner==n: continue
            if n:
                _,j2,r2=primitive_key(n,midpoint,rotations[n]); c2=p[n][j2]
            else:
                r2=r;c2=(arb(0),arb(0))
            v=sub(center,c2)
            equality=frozenset((n,owner))
            for tag,ll in [(k,nl),(j,nh)]:
                if equality<=TIES.get(tag,set()):
                    endpoints+=1
                else:
                    gap=rad-r2+dot(v,ll)
                    assert gap>0,(tag,owner,n,gap)
            # An interior minimum can occur only at normal -v/|v|.
            nv=scale(-1,v)
            if not (det(nl,nv)<0 or det(nv,nh)<0):
                margin=rad-r2-dot(v,v).sqrt()
                assert margin>0,(k,j,owner,n,margin)
                critical+=1
        if j==pieces[active][0]:
            active+=1
    assert active==len(pieces)
    return {'atomic_normal_intervals':len(keys)-1,'exact_contact_endpoint_comparisons':endpoints,
            'interior_minimum_comparisons':critical}


def certify():
    X=[interval(*row) for row in BOX]
    mid=[q((
        F(row[0])+F(row[1]))/2) for row in BOX]
    rad=[q((F(row[1])-F(row[0]))/2) for row in BOX]
    g,H=jets(*X)
    g0,H0=jets(*mid)
    assert H[0][0]>0 and det(*H)>0
    d=det(*H0)
    inv=((H0[1][1]/d,-H0[0][1]/d),(-H0[1][0]/d,H0[0][0]/d))
    C0=tuple(tuple(arb(round(float(v.mid())*10**12))/10**12 for v in row) for row in inv)
    assert det(*C0)>0 or det(*C0)<0
    defect=[[int(i==j)-sum(C0[i][k]*H[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    contraction=max(sum(abs(v).abs_upper() for v in row) for row in defect)
    assert contraction<q('.001')
    for i in range(2):
        deviation=abs(dot(C0[i],g0))+sum(abs(v)*rr for v,rr in zip(defect[i],rad))
        assert deviation<rad[i]
    topo=verify_topology(*X)
    value=area(*mid)
    allowance=sum(abs(dg)*rr for dg,rr in zip(g0,rad))+sum(abs(H[i][j])*rad[i]*rad[j]/2 for i in range(2) for j in range(2))
    value=arb(value,allowance.abs_upper())
    assert q('.83649045546121306')<value<q('.83649045546121308')
    pieces,t,ns,_=chain(*X)
    print(json.dumps({'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'stationary_angle_box':BOX,
        'contraction_norm_upper':str(contraction),'Hessian':[list(map(str,row)) for row in H],
        'Hessian_determinant':str(det(*H)),'placement':list(map(str,t)),
        'support_contact_angles_rad':{str(k):str(ang) for k,ang,*_ in pieces if k in TIES},
        'verified_topology':topo,'area_ball':str(value),
        'strict_area_bracket':['.83649045546121306','.83649045546121308'],
        'scope':'Stationary minimum on the verified two-angle contact manifold, full placement admissible.'},indent=2))

def locate():
    from scipy.optimize import root
    def fun(q):
        grad,hm=jets(*(arb(float(x)) for x in q))
        return list(map(float,grad))
    sol=root(fun,[.512511761647,1.5774790435778],tol=1e-13)
    for x in sol.x:
        print(format(x,'.17g'))
    print('gradient',fun(sol.x),'value',area(*(arb(float(x)) for x in sol.x)))
    print('hessian',jets(*(arb(float(x)) for x in sol.x))[1])
    print('success',sol.success)


if __name__=='__main__':
    parser=ArgumentParser();parser.add_argument('--locate',action='store_true');args=parser.parse_args()
    locate() if args.locate else certify()
