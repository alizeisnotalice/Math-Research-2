"""Enclose the optimum for five regular tests and two smooth width-one bodies inside the unit regular hexagon."""
import json
from fractions import Fraction as Q
import platform

import flint
from flint import arb, ctx
from certify_rotation_inballs import hull

ctx.prec = 128
PI = arb.pi()
GRID = 10**12
N = 512
FOURIER = [
    [(3, Q(1,80)), (5, -Q(19,1200))],
    [(5, Q(1,200)), (15, Q(3,1792))],
]


def exact(q):
    return arb(q.numerator) / q.denominator


def vertex(x, y):
    def grid_point(z):
        k = round(float(z.mid()) * GRID)
        assert abs(z * GRID - k) < 1
        return k
    return grid_point(x), grid_point(y)


def disk_triangle():
    points = []
    for j in range(N):
        s, c = (2*PI*j/N).sin_cos()
        points.append(vertex(c/2, s/2))
    for j in range(3):
        phi = PI/6 + 2*PI*j/3
        s, c = phi.sin_cos()
        x, y = c/arb(3).sqrt(), s/arb(3).sqrt()
        for i in range(N//2+1):
            u, v = (phi+PI-PI/6+2*PI*i/(3*N)).sin_cos()
            points.append(vertex(x+v, y+u))
    return points


def smooth_vertices(terms, phi):
    sin5, cos5 = (5*phi).sin_cos()
    a5 = exact(dict(terms)[5])
    tx, ty = -a5*cos5, a5*sin5
    points = []
    for j in range(N):
        theta = 2*PI*j/N
        sn, cs = theta.sin_cos()
        h, hp = arb(1)/2, arb(0)
        for k, coefficient in terms:
            sine, cosine = (k*(theta-phi)).sin_cos()
            h += exact(coefficient)*cosine
            hp -= k*exact(coefficient)*sine
        points.append(vertex(tx+h*cs-hp*sn, ty+h*sn+hp*cs))
    return points


def regular_vertices(n, phi):
    # Angles are rational multiples of pi, so all arc/vertex choices are exact.
    radius = 1/(2*(PI/(2*n)).cos())
    tx, ty = arb(0), arb(0)
    for j in range(3):
        theta = Q(2*j,3)
        d = abs((theta-phi+Q(1,n)) % Q(2,n) - Q(1,n))
        if d <= Q(1,2*n):
            h = radius*(exact(d)*PI).cos()
        else:
            h = 1-radius*(exact(Q(1,n)-d)*PI).cos()
        sn,cs = (exact(theta)*PI).sin_cos()
        tx += (1-2*h)*cs/3
        ty += (1-2*h)*sn/3
    points = []
    pieces = (N+n-1)//n
    for j in range(n):
        alpha = exact(phi+Q(2*j,n))*PI
        sn,cs = alpha.sin_cos()
        x,y = tx+radius*cs, ty+radius*sn
        for i in range(pieces+1):
            normal = alpha+PI*(1-exact(Q(1,2*n))+exact(Q(i,n*pieces)))
            sn,cs = normal.sin_cos()
            points.append(vertex(x+cs,y+sn))
    return points


def bounds(poly, delta):
    twice = sum(x*b-y*a for (x,y),(a,b) in zip(poly,poly[1:]+poly[:1]))
    area = Q(twice,2*GRID**2)
    eta = Q(2,GRID)
    # P is inside the regular hexagon enlarged by eta; its perimeter is <4.
    return area/(1+2*eta)**2, area+4*delta+Q(22,7)*delta**2


def data():
    base = hull(disk_triangle())
    groups = []
    for n,count in ((5,6),(7,6),(9,2)):
        denom = 6*n if n != 9 else 2*n
        groups.append([hull(regular_vertices(n,Q(2*j+1,denom))) for j in range(count)])
    groups += [[hull(smooth_vertices(terms,(2*j+1)*PI/(2*k))) for j in range(6)]
               for k,terms in zip((3,15),FOURIER)]
    curvature = max(Q(1,2)+sum((k+1)*(k*k-1)*abs(a) for k,a in terms)
                    for terms in FOURIER)
    delta = Q(2,GRID)+max(Q(1),curvature)*(2*Q(22,7)/N)**2/8
    return base,groups,delta


def certify():
    global N
    base,groups,delta = data()
    incumbent = None
    candidate_intervals = []
    def descend(poly, choices):
        nonlocal incumbent
        if len(choices) == len(groups):
            lower,upper = bounds(poly,delta)
            candidate_intervals.append((lower,upper,choices))
            if incumbent is None or upper < incumbent:
                incumbent = upper
            return
        for i,points in enumerate(groups[len(choices)]):
            descend(hull(poly+points),choices+(i,))
    descend(base,())
    coarse_lower = min(row[0] for row in candidate_intervals)
    promising = [row for row in candidate_intervals if row[0]<=incumbent]
    coarse_n = N
    N = 4096
    base,groups,delta = data()
    refined = []
    for _,_,indices in promising:
        poly = hull(base+[pt for choices,i in zip(groups,indices) for pt in choices[i]])
        lo,hi = bounds(poly,delta)
        refined.append((lo,hi,indices))
    omitted = [row[0] for row in candidate_intervals if row[0]>incumbent]
    return {"python":platform.python_version(),"python_flint":flint.__version__,
            "grid":GRID,"coarse_samples_per_turn":coarse_n,"fine_samples_per_turn":N,
            "precision":ctx.prec,"configurations":len(candidate_intervals),
            "coarse_lower":str(coarse_lower),"coarse_upper":str(incumbent),
            "omitted_lower":str(min(omitted)) if omitted else None,
            "lower":str(min([r[0] for r in refined]+omitted)),
            "upper":str(min(r[1] for r in refined)),
            "refined":[[list(ids),str(a),str(b)] for a,b,ids in refined]}


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2))
