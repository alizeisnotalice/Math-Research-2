"""Arb bounds for the least parallel padding of the four-body minimizing hull."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform

import flint
from flint import arb, ctx

from certify_four_body_basin import CENTERS
from certify_four_body_contact_stationary import (
    BOX, add, chain, det, dot, geometry, interval, normal, pi, q, r, regular, sub, verify_topology,
)
from certify_four_body_translation_basin import first
from four_contact_taylor import Model

TURN = 2*pi/9
PAD_LO = F('.00394138015')
PAD_HI = F('.00394138020')
MOVE = ('.35681038', '.003790506253', '-.001905119134')
SEED_DENOM = 10**12
AREA_THRESHOLD = F('.836490455462')
AREA_LOWER = F('.83649045546121306')
STABLE_PADDING_LOWER = F('.00392411')
HOST_SOURCE = Path(__file__).with_name('certify_four_body_contact_stationary.py')


def axis(seed):
    sign, key = seed
    assert type(sign) is int and sign in (-1,1)
    assert type(key) is int and abs(key) <= SEED_DENOM
    s=F(key,SEED_DENOM)
    return (sign*(1-s*s)/(1+s*s),sign*2*s/(1+s*s)),2*q(s).atan()+(0 if sign==1 else pi)


def balance(directions):
    a,b,c=directions
    w=det(b,c),det(c,a),det(a,b)
    assert all(x>0 for x in w) or all(x<0 for x in w)
    ws=tuple(x/sum(w) for x in w)
    assert min(ws)>0 and sum(ws)==1
    assert all(sum(l*d[k] for l,d in zip(ws,directions))==0 for k in (0,1))
    return ws


def features(n, angle, rotation):
    """Potential support vertex or arc centre over a small real angular ball."""
    gap=pi/(2*n)
    step=2*pi/n
    out=[]
    for radius in (0,1):
        center=rotation+radius*pi
        index=round(float(((angle-center)/step).mid()))
        for k in (index-1,index,index+1):
            delta=(angle-center)-k*step
            if not abs(delta)>gap:
                out.append((k%n,radius,abs(delta)<gap))
        # This excludes any eligible index outside the three tested above.
        assert abs((angle-center)-index*step)<step
    assert out
    return out


class Hull:
    def __init__(self):
        root = [interval(*b) for b in BOX]
        verify_topology(*root)
        self.pieces, self.position, _, self.vertices = chain(*root)
        self.rot={3:arb(0),5:self.position[2],7:self.position[5]}
        self.seed_cache={}

    def support(self,seed):
        if seed in self.seed_cache:
            return self.seed_cache[seed]
        u,angle=axis(seed)
        v=tuple(map(q,u))
        upper=r
        for n in (3,5,7):
            for k,radius,_ in features(n,angle,self.rot[n]):
                upper=max(upper,(dot(v,self.vertices[n][k])+radius).upper())
        self.seed_cache[seed]=u,v,angle,upper
        return self.seed_cache[seed]


def material(u, direction, theta):
    """One support material point at the angular centre, checked in the closed body."""
    options=features(9,direction,theta)
    k,arc,strict=options[0]
    if strict:
        x=regular(9,k)
        R=normal(theta)
        p=(R[0]*x[0]-R[1]*x[1],R[1]*x[0]+R[0]*x[1])
        return add(p,u) if arc else p
    # The normal may be exactly on a switch. Its nearest vertex is valid material.
    j=round(float(((direction-theta)/TURN).mid()))
    x=r/(pi/18).cos()
    a=theta+j*TURN
    return x*a.cos(),x*a.sin()


def lower_cell(host, depth, index, seeds):
    s=F(2*index+1,2**(depth+1))
    theta=q(s)*TURN
    width=TURN/2**(depth+1)
    data=[host.support(tuple(seed)) for seed in seeds]
    weights=balance([v[0] for v in data])
    coeff=[arb(0) for _ in range(3)]
    for lam,(_,u,angle,h) in zip(weights,data):
        z=material(u,angle,theta)
        for k,value in enumerate((dot(z,u),det(z,u),h)):
            coeff[k]+=q(lam)*value
    A,B,C=coeff
    assert A>0 and width>0 and width<pi/2
    return (A*width.cos()-abs(B)*width.sin()-C).lower()


def upper_bound(host):
    theta,tx,ty=map(q,MOVE)
    t=(tx,ty)
    cuts=[arb(0),*(a for _,a,_,_,_ in host.pieces)]
    for k,sgn in product(range(9),(-1,1)):
        alpha=theta+k*TURN+sgn*TURN/4
        if alpha<0:
            alpha+=2*pi
        if alpha>2*pi:
            alpha-=2*pi
        assert 0<alpha<2*pi
        cuts.append(alpha)
    cuts.sort(key=lambda a:float(a.mid()))
    assert all(a<b for a,b in zip(cuts,cuts[1:]))
    local = []
    host_i=0
    length=0
    for l,h in zip(cuts,cuts[1:]):
        m=(l+h)/2
        while m>host.pieces[host_i][1]:
            host_i+=1
        prior=arb(0) if host_i==0 else host.pieces[host_i-1][1]
        assert prior<m<host.pieces[host_i][1] and h-l<pi
        _,_,_,host_radius,c=host.pieces[host_i]
        b9 = features(9,m,theta)
        assert len(b9)==1 and b9[0][2]
        k,radius,_=b9[0]
        p=tuple(r/(pi/18).cos()*v for v in normal(theta+k*TURN))
        d=sub(add(t,p),c)
        shift=radius-host_radius
        ul,uh=normal(l),normal(h)
        local.extend([(shift+dot(d,ul)).upper(),(shift+dot(d,uh)).upper()])
        # If d is not outside the wedge, admit its direction as a possible maximum.
        if not (det(ul,d)<0 or det(d,uh)<0):
            local.append((shift+(dot(d,d)).sqrt()).upper())
        length+=1
    return max(local),length


def replay(certificate):
    host=Hull()
    assert certificate['strict_padding_bracket']==[str(PAD_LO),str(PAD_HI)]
    assert certificate['stationary_source_sha256']==sha256(HOST_SOURCE.read_bytes()).hexdigest()
    assert certificate['rotation_and_translation']==list(MOVE)
    cells=certificate['cover']
    endpoints=[]
    lowers=[]
    for depth,index,seeds in cells:
        assert type(depth) is int and 0<=depth<=50
        assert type(index) is int and 0<=index<2**depth
        assert len(seeds)==3
        endpoints.append((F(index,2**depth),F(index+1,2**depth)))
        lb=lower_cell(host,depth,index,seeds)
        assert lb>q(PAD_LO), (depth,index,lb)
        lowers.append(lb)
    endpoints.sort()
    assert endpoints[0][0]==0 and endpoints[-1][1]==1
    assert all(a[1]==b[0] for a,b in zip(endpoints,endpoints[1:]))
    upper,interval_count=upper_bound(host)
    assert upper<q(PAD_HI),upper
    return {'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'stable_padding':stable_padding(),
            'cells':len(cells),'deepest':max(row[0] for row in cells),'host_directions':len(host.seed_cache),
            'lower_floor':str(min(lowers)),'candidate_upper':str(upper),'upper_atomic_intervals':interval_count}


def stable_padding():
    """Quantitative inheritance from the low-area translation basin."""
    budget = AREA_THRESHOLD-AREA_LOWER
    angle_offset = F('.000003968')
    translation_offset = F('.000000004373')
    assert angle_offset**2 > 20*budget
    preliminary = F('.00000001431')
    assert preliminary*F('.000055')>budget
    assert translation_offset*(F('.00018')-F(5,4)*preliminary)>budget

    Model.degree=8
    R=q('.00001')
    h=arb(0,q('.4'))
    centers=list(map(q,CENTERS))
    assert all((interval(*row)-c).abs_upper()+q(angle_offset)<q('.000004')
               for c,row in zip(centers,BOX))
    alpha,beta=(Model({(0,0):c,k:R}) for c,k in zip(centers,[(1,0),(0,1)]))
    position=geometry(alpha,beta)[0]
    derivative_squares=[sum(first(e,j,h,h,R).abs_upper()**2 for j in (0,1))
                        for e in position]
    bounds=[]
    for n,k,i in [(3,0,None),(5,3,2),(7,6,5)]:
        z=(derivative_squares[k]+derivative_squares[k+1]).sqrt()
        if i is not None:
            z+=r/(pi/(2*n)).cos()*derivative_squares[i].sqrt()
        bounds.append(z.abs_upper())
    L=F('4.351')
    assert all(z<q(L) for z in bounds)
    hausdorff_upper=L*angle_offset+translation_offset
    assert STABLE_PADDING_LOWER<PAD_LO-hausdorff_upper
    return {'area_threshold':str(AREA_THRESHOLD),'tangent_offset_upper':str(angle_offset),
            'translation_offset_l1_upper':str(translation_offset),'body_lipschitz_upper':list(map(str,bounds)),
            'hausdorff_upper':str(hausdorff_upper),'stable_padding_lower':str(STABLE_PADDING_LOWER)}


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('certificate',type=Path)
    args=parser.parse_args()
    raw=args.certificate.read_bytes()
    result=replay(json.loads(raw))
    result['certificate_sha256']=sha256(raw).hexdigest()
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
