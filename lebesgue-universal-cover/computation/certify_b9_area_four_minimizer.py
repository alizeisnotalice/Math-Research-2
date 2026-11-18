"""Integer and Arb verification of B9's area cost on the four-body hull."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import platform

import flint
from flint import arb, ctx

from certify_four_body_basin import CENTERS, RAW_CENTER
from certify_four_body_contact_stationary import BOX, geometry, interval, pi, r
from certify_four_body_translation_basin import first
from four_contact_taylor import Model
from certify_five_harmonic import disk_points, marked_hull, sqrt_upper
from certify_nine import S, PI_UPPER, append_point, ar, arc_body, arithmetic, hull

MOVE = ('.349106707', '-.000764981', '-.000579172')
UPPER = F('.83730940')
T_LIMIT = F(1,4)
BODY_DATA = (
    (3, None, 0, F('.42'), F('.29')),
    (5, 2, 3, F('.47'), F('.163')),
    (7, 5, 6, F('.48'), F('.115')),
    (9, 8, 9, F('.49'), F('.089')),
)
POSE_EPS = F(1,10**14)
CASES = {
    'window':(192,192,F('.83724')),
    'fixed':(640,640,F('.8373092')),
    'sublevel':(384,384,F('.837297')),
}


def validate():
    assert pi<ar(PI_UPPER)
    for n,_,_,inner,L in BODY_DATA:
        gap=pi/(2*n)
        assert 1-1/(2*gap.cos())>ar(inner)
        assert gap.tan()/2<ar(L)
    actual=geometry(*(interval(*v) for v in BOX))[0]
    assert all(abs(a-ar(c))<ar(POSE_EPS) for a,c in zip(actual,RAW_CENTER))


def sublevel_errors():
    """Support perturbations from the radius and displacement of the near-minimizer."""
    budget = F('.836490455462')-F('.83649045546121306')
    angle = F('.000003968')
    shift = F('.000000004373')
    preliminary = F('.00000001431')
    assert angle**2 > 20*budget
    assert preliminary*F('.000055')>budget
    assert shift*(F('.00018')-F(5,4)*preliminary)>budget

    Model.degree=8
    R=ar(F('.00001'))
    h=arb(0,ar(F('.4')))
    C=list(map(ar,map(F,CENTERS)))
    assert all((interval(*x)-c).abs_upper()+ar(angle)<ar(F('.000004')) for c,x in zip(C,BOX))
    models=[Model({(0,0):c,k:R}) for c,k in zip(C,[(1,0),(0,1)])]
    t=geometry(*models)[0]
    J=[[first(v,j,h,h,R) for j in (0,1)] for v in t]
    lipschitz=list(map(F,('.402','1.33','1.5')))
    maxima=[]
    for (n,i,k,_,_),L in zip(BODY_DATA,lipschitz):
        radius=r/(pi/(2*n)).cos()
        rotation=arb(0) if i is None else t[i].evaluate(h,h)
        upper=arb(0)
        for j in range(n):
            phi=rotation+2*pi*j/n
            for v in range(128):
                tau=arb(pi*(2*v+1-128)/(2*n*128),pi/(2*n*128))
                U=phi+tau
                co,si=U.cos(),U.sin()
                rotation_derivative=radius*tau.sin()
                row=[J[k][l]*co+J[k+1][l]*si+(rotation_derivative*J[i][l] if i is not None else 0)
                     for l in (0,1)]
                norm=sum(a.abs_upper()**2 for a in row).sqrt().upper()
                upper=max(upper,norm)
        assert upper<ar(L),(n,upper,L)
        maxima.append(str(upper))
    errors=[L*angle+shift+3*POSE_EPS for L in lipschitz]
    return errors,{'body_support_lipschitz':maxima,'support_from_reference':list(map(str,errors))}


def common_host(k, window=None, support_errors=None):
    points=[]
    disk_points(points,8*k,pi,F(1,2)-F(2,S))
    assert all(x*x+y*y<=(S//2)**2 for x,y in points)
    half=[POSE_EPS]*8
    if window is not None:
        d=window['denominator']
        assert d==2**44
        assert all(len(window[x])==8 for x in ('low','high'))
        assert all(type(v) is int for side in ('low','high') for v in window[side])
        assert all(F(l,d)<c<F(h,d) for l,h,c in zip(window['low'],window['high'],RAW_CENTER))
        half=[max(c-F(l,d),F(h,d)-c) for l,h,c in zip(window['low'],window['high'],RAW_CENTER)]
    shrink=[]
    for n,j,k0,inner,L in BODY_DATA[:3]:
        phi=RAW_CENTER[j] if j is not None else F(0)
        w=half[j] if j is not None else F(0)
        a,b=half[k0:k0+2]
        e=(sqrt_upper(a*a+b*b)+L*w if support_errors is None else support_errors[len(shrink)])+F(4,S)
        factor=1-e/inner
        assert 0<factor<1
        body=[]
        arc_body(body,k,pi,n,ar(phi),arb(0),arb(0))
        t=tuple(c*S for c in RAW_CENTER[k0:k0+2])
        points.extend((round(factor*x+t[0]),round(factor*y+t[1])) for x,y in body)
        shrink.append(str(factor))
    return hull(points),shrink


class AngleCell:
    def __init__(self, old, material, depth, index):
        center=2*pi*ar(F(2*index+1,9*2**(depth+1)))
        w=2*PI_UPPER*F(1,9*2**(depth+1))
        c,s=center.cos(),center.sin()
        _,_,_,inner,L=BODY_DATA[-1]
        factor=1-(L*w+F(5,S))/inner
        assert 0<factor<1
        new=[]
        for x,y in material:
            x,y=ar(F(x,S)),ar(F(y,S))
            append_point(new,c*x-s*y,s*x+c*y)
        self.p=[(x,y,0) for x,y in old]
        self.q=hull((round(factor*x),round(factor*y)) for x,y in new)

    def plane(self, site):
        tx,ty=site
        H=marked_hull(self.p+[(x+tx,y+ty,1) for x,y in self.q])
        v=gx=gy=0
        for (x,y,b),(u,w,c) in zip(H,H[1:]+H[:1],strict=True):
            v+=x*w-y*u
            gx+=b*w-c*y
            gy+=c*x-b*u
        return (v-gx*tx-gy*ty,gx,gy)


def affine_floor(planes, weights):
    assert all(type(x) is int and x>0 for x in weights)
    assert len(weights)==len(planes)
    v,x,y=[sum(w*p[j] for p,w in zip(planes,weights,strict=True)) for j in range(3)]
    return F(v-T_LIMIT*S*(abs(x)+abs(y)),2*S*S*sum(weights))


def feasible():
    k=4096
    points=[]
    disk_points(points,8*k,pi,F(1,2))
    actual=geometry(*(interval(*v) for v in BOX))[0]
    pose=list(actual)+list(map(ar,MOVE))
    for n,angle,k0,_,_ in BODY_DATA:
        arc_body(points,k,pi,n,pose[angle] if angle is not None else arb(0),*pose[k0:k0+2])
    result=arithmetic(hull(points),k)
    eps=PI_UPPER**2/(72*k*k)+F(2,S)
    polygon_area=F(result['grid_polygon_area'])
    upper=polygon_area+F(result['hull_perimeter_upper'])*eps+PI_UPPER*eps**2
    assert upper<UPPER, float(upper)
    return {'pose':MOVE,'arc_subdivisions':k,'strict_upper':str(upper),
            'strict_upper_display':float(upper),'polygon_area':str(polygon_area),
            'hull_perimeter_upper':result['hull_perimeter_upper'],'arc_support_allowance':str(eps)}


def outside():
    z=ar(BODY_DATA[-1][3]+T_LIMIT)
    far=(pi-(1/(2*z)).acos())/4+(z*z-arb(1)/4).sqrt()/2
    assert far>ar(UPPER)
    return far


def verify(certificate, window):
    validate()
    far=outside()
    report={}
    assert set(certificate['cases'])==set(CASES)
    errors,summary=sublevel_errors()
    for name,(k,k9,target) in CASES.items():
        old,shrink=common_host(k, window if name=='window' else None, errors if name=='sublevel' else None)
        new=[]
        arc_body(new,k9,pi,9,arb(0),arb(0),arb(0))
        new=hull(new)
        bounds=[]
        floors=[]
        for depth,index,rows in certificate['cases'][name]:
            assert type(depth) is int and 0<=depth<=35
            assert type(index) is int and 0<=index<2**depth
            assert rows and all(len(a)==3 and all(type(v) is int for v in a) for a in rows)
            assert all(max(abs(x),abs(y))<=S*T_LIMIT for x,y,_ in rows)
            bounds.append((F(index,2**depth),F(index+1,2**depth)))
            cell=AngleCell(old,new,depth,index)
            f=affine_floor([cell.plane((x,y)) for x,y,_ in rows],[w for _,_,w in rows])
            assert f>target, (name,depth,index,float(f))
            floors.append(f)
        bounds.sort()
        assert bounds[0][0]==0 and bounds[-1][1]==1 and all(a[1]==b[0] for a,b in zip(bounds,bounds[1:]))
        floor=min(floors)
        assert far>ar(floor)
        report[name]={'lower':str(floor),'lower_display':float(floor),'cells':len(floors),
                      'deepest':max(row[0] for row in certificate['cases'][name]),
                      'host_subdivisions':k,'b9_subdivisions':k9,'body_shrink':shrink,'host_vertices':len(old),
                      'threshold':str(target)}
    return {'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
            'far_disk_point_lower':str(far.lower()),'sublevel':summary,'feasible':feasible(),'cases':report}


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('certificate',type=Path)
    parser.add_argument('dyadic_receipt',type=Path)
    args=parser.parse_args()
    raw=args.certificate.read_bytes()
    receipt=args.dyadic_receipt.read_bytes()
    assert sha256(receipt).hexdigest()=='4086b16e625cc2871c7e8c34b3802681857fd0ba21ddfafc80b2cb0aa1ce51c5'
    result=verify(json.loads(raw),json.loads(receipt)['exact_dyadic_window'])
    result['source_sha256']={p.name:sha256(p.read_bytes()).hexdigest() for p in
        (Path(__file__),args.certificate,args.dyadic_receipt)}
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
