"""Fixed finite positive-part area duals for C4 and a regular pentagon."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import sys
import flint
from flint import arb, ctx

ctx.prec = 192
pi = arb.pi()
half = arb(1)/2
root = arb(3).sqrt()
wall_length = 1/(2*root)
r = half/(pi/10).cos()

def rad(x):
    return pi*x/180

def direction(theta):
    return theta.cos(),theta.sin()

def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]

def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]

def scale(s,v):
    return s*v[0],s*v[1]

def bounds(x):
    def rational(y):
        m,k=y.man_exp()
        return str(Fraction(int(m))*Fraction(2)**int(k))
    return [rational(x.lower()),rational(x.upper())]

# Each row contains (pentagon vertex, total curvature measure, first moment).
def arc(lo,hi,j):
    return j,(hi-lo)/2,((hi.sin()-lo.sin())/2,(lo.cos()-hi.cos())/2)

def wall(a,j):
    return j,wall_length,scale(wall_length,direction(rad(a)))

def row(a,b,j):
    if a==b:
        assert a in (0,60,120,180,240,300)
        return wall(a,j)
    assert any(s<=a<b<=s+60 for s in (60,180,300))
    return arc(rad(a),rad(b),j)

def balanced(full,partial):
    x,y=partial
    target=tuple(-sum(v[k] for _,_,v in full) for k in (0,1))
    determinant=cross(x[2],y[2])
    weights=(cross(target,y[2])/determinant,cross(x[2],target)/determinant)
    assert all(0<w<1 for w in weights)
    blocks=full+[(j,w*m,scale(w,u)) for w,(j,m,u) in zip(weights,partial)]
    return blocks,weights

def support_term(blocks,theta):
    return sum(r*dot(direction(theta+rad(72*j)),v) for j,_,v in blocks)

def replay():
    p=2*(r*rad(6).cos()-half)/root
    top=r+p
    q3=(r*rad(234).cos(),r*rad(234).sin()+p)
    a=(half/top).acos()
    b=(half/dot(q3,q3).sqrt()).acos()
    inner=pi/8+root/4
    exact=inner+((a.tan()-a)+2*(b.tan()-b))/4
    assert 24*pi/180<a<26*pi/180
    assert 7*pi/180<b<8*pi/180
    assert r*rad(36).cos()+p<half
    assert top*rad(30).cos()<half
    assert cross(q3,direction(rad(240)))>0
    assert 0<r*rad(162).sin()+p
    # Arc unit vectors and the corresponding maximizing vertices.
    exact_rows=[arc(pi/2-a,pi/2+a,1),arc(rad(240)-2*b,rad(240),3),
                arc(rad(300),rad(300)+2*b,4)]
    data=[
        ((6,11),[row(0,0,0),row(60,105,1),row(202,240,3)],
                   [row(240,240,3),row(359,360,0)]),
        ((11,15),[row(240,240,3),row(63,112,1),row(214,240,3)],
                    [row(0,0,0),row(62,63,1)]),
        ((15,17),[row(64,115,1),row(222,240,3),row(301,312,4)],
                    [row(0,0,0),row(240,240,3)]),
        ((17,18),exact_rows,[wall(240,3),wall(0,0)]),
    ]
    certificates=[]
    for i,((lo,hi),full,partial) in enumerate(data):
        blocks,weights=balanced(full,partial)
        total_mass=sum(m for _,m,_ in blocks)
        def lower(theta):
            return inner+(support_term(blocks,theta)-half*total_mass)/2
        lob=lower(rad(lo))-exact
        hib=lower(rad(hi))-exact
        concavity=support_term(blocks,rad(lo).union(rad(hi)))/2
        assert lob>0 and concavity>0
        if i<3:
            assert hib>0
        else:
            assert hib.contains(0)  # Equality itself is proved by exact contacts.
        certificates.append(dict(interval=[lo,hi],weights=[bounds(w) for w in weights],
                                 concavity_margin=bounds(concavity),lower_end_surplus=bounds(lob),
                                 upper_end_surplus=bounds(hib) if i<3 else 'exact equality'))
    print('Minimum area',exact)
    print('Increment',exact-inner)
    print('min endpoint margin',min(float(Fraction(d['lower_end_surplus'][0])) for d in certificates))
    return dict(exact_area=bounds(exact),area_increment=bounds(exact-inner),initial_area=bounds(inner),
                angle_of_pentagon_degrees=18,translation_x=0,translation_y=bounds(p),
                exposed_half_angles_radians=[bounds(v) for v in (a,b)],certificates=certificates,
                python=sys.version,python_flint=flint.__version__,precision_bits=ctx.prec,
                program_sha256=sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    Path(sys.argv[1]).write_text(json.dumps(replay(),indent=2)+'\n')
