"""Contact-coordinate basin with fixed angles and translated material fans."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb, arb_mat, ctx

from certify_four_body_basin import CENTERS, RAW_CENTER, area_model, fan_parts, mapped_data, topology
from certify_four_body_contact_stationary import BOX, chain, det, interval, normal, q, sub
from four_contact_taylor import Model, powers

INDICES = {3: 0, 5: 2, 7: 4}
TRANSLATIONS = (0, 1, 3, 4, 6, 7)
ROTATIONS = (2, 5)
K = F('.0004')
MU = F('.1')
DEGREE = 16


def first(model, variable, x, y, scale):
    value = arb(0)
    xp,yp=powers(x,Model.degree),powers(y,Model.degree)
    for (i, j), a in model.poly.items():
        power = (i, j)[variable]
        if power:
            value += a*power*xp[i-(variable==0)]*yp[j-(variable==1)]
    margin = 1 - max(x.abs_upper(), y.abs_upper())
    assert margin > 0
    extra = model.error/margin
    return arb(value, extra.abs_upper())/scale


def precondition(matrix):
    m = arb_mat([[entry.evaluate(arb(0), arb(0)) for entry in row] for row in matrix])
    inv = m.inv()
    return arb_mat([[q(F(round(float(inv[i,j].mid())*10**12), 10**12))
                     for j in range(m.ncols())] for i in range(m.nrows())])


def material_fans(data):
    arc_grad = [Model() for _ in INDICES for _ in (0,1)]
    start = Model()
    for _, stop, owner, radius, _ in data[0]:
        if owner and radius:
            diff = sub(normal(stop), normal(start))
            j = INDICES[owner]
            arc_grad[j] += radius*diff[1]/2
            arc_grad[j+1] -= radius*diff[0]/2
        start = stop
    for _, links, _, _ in fan_parts(data, correlated=True):
        g = arc_grad.copy()
        h = [[0]*6 for _ in range(6)]
        edges = []
        for (owner_p,P), (owner_q,Q) in links:
            if owner_p:
                p = INDICES[owner_p]
                g[p] += Q[1]/2
                g[p+1] -= Q[0]/2
            if owner_q:
                j = INDICES[owner_q]
                g[j] -= P[1]/2
                g[j+1] += P[0]/2
            if owner_p and owner_q:
                assert p != j
                h[p][j+1] += 1
                h[j+1][p] += 1
                h[p+1][j] -= 1
                h[j][p+1] -= 1
            edges.append((bool(owner_p),bool(owner_q), P,Q,det(P,Q)))
        h_norm = F(max(sum(abs(x) for x in row) for row in h),2)
        yield g, h_norm, edges


def process_fans(data):
    fan_data = list(material_fans(data))
    G = [[Model(1)]*7, *[[g[i] for g,_,_ in fan_data] for i in range(6)]]
    inv = precondition(G)
    assert abs(inv.det()) > 0
    return fan_data, inv, [[sum(inv[i,k]*G[k][j] for k in range(7))
                            for j in range(7)] for i in range(7)]


def orientation_floor(fan_data, x, y, D):
    def size(P):
        return sum(z.evaluate(x,y).abs_upper() for z in P)
    margins=[]
    for _,_,edges in fan_data:
        for p,j,P,Q,C in edges:
            lower=C.evaluate(x,y).lower() - D*(p*size(Q)+j*size(P)) - (2*p*j)*D**2
            margins.append(lower.lower())
    return min(margins)


def cell(data, value, fans, D, scale, x, y, K):
    H = [[value.second(i,j,x,y,scale) for j in range(2)] for i in range(2)]
    for i in range(2):
        H[i][i] -= q(MU)
    assert H[0][0]>0 and det(*H)>0, ('curvature', H[0][0], det(*H))
    topology(mapped_data(data,x,y))
    fan_data, B, pre_G = fans
    M = arb_mat([[entry.evaluate(x,y) for entry in row] for row in pre_G]).inv()*B
    floor = min((M[i,0]-K*sum(abs(M[i,j]) for j in range(1,7))).lower() for i in range(7))
    assert floor>0, ('simplex', floor)
    orientation = orientation_floor(fan_data,x,y,D)
    assert orientation>0, ('orientation',orientation)
    return (H[0][0].lower(), det(*H).lower(), floor, orientation)


def localize_angles(placements, scale, target, rect):
    """Find one contact anchor for every raw-angle pair by a contraction."""
    m = [placements[k] for k in ROTATIONS]
    J0 = [[first(e,j,arb(0),arb(0),scale) for j in range(2)] for e in m]
    J = arb_mat(J0).inv()
    B = [[q(F(round(float(J[i,j].mid())*10**12),10**12)) for j in range(2)] for i in range(2)]
    BD = [sum(B[i][j]*m[j] for j in range(2)) for i in range(2)]
    assert abs(det(*B)) > 0
    x=arb(0,rect)
    M = [[int(i==j)-first(BD[i],j,x,x,scale) for j in range(2)] for i in range(2)]
    lipschitz = max(sum(a.abs_upper() for a in row) for row in M)
    assert lipschitz<1, ('inverse_derivative',lipschitz)
    V = [sum(B[i][j]*(target[j]-m[j].evaluate(arb(0),arb(0))) for j in range(2))
         for i in range(2)]
    assert all(v.abs_upper()+lipschitz*rect*scale<rect*scale for v in V), ('inverse_angles',V,lipschitz)
    u = [arb(v,sum(a.abs_upper() for a in row)*rect*scale) for v,row in zip(V,M)]
    for _ in range(5):
        u = [V[i]+sum(M[i][j]*u[j] for j in range(2)) for i in range(2)]
    return [a/scale for a in u],lipschitz


def source_ends(source, halfwidths):
    if source:
        raw = json.loads(source.read_bytes())
        return [[F(a,2**44) for a in raw[k]] for k in ('L','H')]
    return [[p+s*w for p,w in zip(RAW_CENTER,halfwidths)] for s in (-1,1)]


def main():
    parser = ArgumentParser()
    parser.add_argument('--radius', type=F, default=F('.0002'))
    parser.add_argument('--transverse', type=F, default=F('.0002'))
    parser.add_argument('--gradient-cube', type=F, default=K)
    parser.add_argument('--degree', type=int, default=DEGREE)
    parser.add_argument('--depth', type=int, default=14)
    parser.add_argument('--raw', type=F, nargs=8)
    parser.add_argument('--endpoints', type=Path)
    args = parser.parse_args()
    Model.degree = args.degree
    A,B = map(q,CENTERS)
    R = q(args.radius)
    D = q(args.transverse)
    alpha=Model({(0,0):A,(1,0):R}); beta=Model({(0,0):B,(0,1):R})
    data = chain(alpha,beta)
    half=F(1,2)
    assert all(abs((interval(*s)-c)/R)<q(half) for s,c in zip(BOX,(A,B)))
    results={}; counts=Counter()
    if args.raw or args.endpoints:
        lo,hi = source_ends(args.endpoints,args.raw)
        raw = [interval(a,b) for a,b in zip(lo,hi)]
        uv,lip = localize_angles(data[1],R,[raw[i] for i in ROTATIONS],q(half))
        assert all(a.abs_upper()<q(half) for a in uv)
        diff = [(raw[i]-data[1][i].evaluate(*uv)).abs_upper() for i in TRANSLATIONS]
        assert all(v<D for v in diff), ('raw_translations',diff)
        results = {'raw_lo':list(map(str,lo)), 'raw_hi':list(map(str,hi)),
                   'raw_contraction':str(lip), 'normalized_tangent':list(map(str,uv)),
                   'translation_offsets':list(map(str,diff))}
    value = area_model(data[0])
    fan_data, _, _ = fans = process_fans(data)
    remainder = max(h_norm for _,h_norm,_ in fan_data)*args.transverse/2
    assert remainder < args.gradient_cube
    worst = None
    def visit(lo,hi,depth):
        nonlocal worst
        counts['visited']+=1
        x,y = [interval(l,h) for l,h in zip(lo,hi)]
        try:
            margins = cell(data,value,fans,D,R,x,y,q(args.gradient_cube))
        except (AssertionError,ZeroDivisionError,ValueError) as e:
            if depth >= args.depth:
                raise
            counts[e.args[0][0] if e.args and isinstance(e.args[0],tuple) else str(e)[:80]]+=1
            j=depth%2
            m=(lo[j]+hi[j])/2
            visit(lo,hi[:j]+(m,)+hi[j+1:],depth+1)
            visit(lo[:j]+(m,)+lo[j+1:],hi,depth+1)
            return
        counts['terminal']+=1
        worst=margins if worst is None else tuple(min(a,b) for a,b in zip(worst,margins))
    visit((-half,-half),(half,half),0)
    assert counts['visited']==2*counts['terminal']-1
    source=Path(__file__)
    results.update({'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
        'script_sha256':sha256(source.read_bytes()).hexdigest(),'model_degree':args.degree,
        'area_model_remainder':str(value.error),'contact_center':CENTERS, 'contact_halfwidth':str(args.radius/2),
        'translation_halfwidth':str(args.transverse),'gradient_cube':str(args.gradient_cube),
        'curvature_lower':str(MU),'linear_remainder':str(remainder),'fan_margin':str(args.gradient_cube-remainder),
        'counts':counts,'min_margins':list(map(str,worst))})
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    main()
