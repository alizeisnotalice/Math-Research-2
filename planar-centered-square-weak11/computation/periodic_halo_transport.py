#!/usr/bin/env python3
"""Full periodic halo with a verified label-preserving transport tail bound.

Derived from periodic_halo.py SHA256
b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33.
"""
import argparse, bisect, hashlib, json, math, platform, sys, time
from fractions import Fraction as F
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def union_area(rects):
    if not rects:
        return 0
    ys = sorted({v for r in rects for v in (r[2], r[3])})
    yi = {y: i for i, y in enumerate(ys)}
    events = sorted((x, sign, yi[y0], yi[y1]) for x0, x1, y0, y1 in rects for x, sign in ((x0, 1), (x1, -1)))
    n = len(ys) - 1
    cover, length = [0] * (4*n+4), [0] * (4*n+4)
    def update(node, lo, hi, a, b, delta):
        if a <= lo and hi <= b:
            cover[node] += delta
        else:
            mid = (lo + hi)//2
            if a < mid: update(node*2, lo, mid, a, b, delta)
            if b > mid: update(node*2+1, mid, hi, a, b, delta)
        if cover[node]: length[node] = ys[hi] - ys[lo]
        elif hi-lo == 1: length[node] = 0
        else: length[node] = length[node*2] + length[node*2+1]
    area = 0
    previous = events[0][0]
    for x, sign, a, b in events:
        area += (x-previous)*length[1]
        update(1, 0, n, a, b, sign)
        previous = x
    return area


def calculate_transport(data, base, copy_factors, delta, bits=40, R=F(5), labels=None, rectangles_out=None):
    tail, delta = F(R), F(delta)
    if len(copy_factors) != 2 or any(type(v) is not int or v <= 0 for v in copy_factors):
        raise ValueError('copy factors must be positive integers')
    nxcopy, nycopy = copy_factors
    bx, by, blam = (F(base[k]) for k in ('Px','Py','lam'))
    px, py, lam = (F(data[k]) for k in ('Px','Py','lam'))
    if min(bx,by,blam,tail) <= 0 or delta < 0 or (px,py,lam) != (nxcopy*bx,nycopy*by,blam):
        raise ValueError('periods, threshold, radius or delta invalid')
    source = [(F(x)%bx,F(y)%by,F(w)) for x,y,w in base['pts']]
    actual = [(F(x)%px,F(y)%py,F(w)) for x,y,w in data['pts']]
    if not source or any(w < 0 for _,_,w in source):
        raise ValueError('base requires nonnegative masses and positive total')
    total = sum((w for _,_,w in source),F(0))
    if total <= 0 or len(actual) != nxcopy*nycopy*len(source):
        raise ValueError('base mass or number of matched rows invalid')
    expected_labels = {(a,b,i) for a in range(nxcopy) for b in range(nycopy) for i in range(len(source))}
    if labels is None:
        labels = [(a,b,i) for a in range(nxcopy) for b in range(nycopy) for i in range(len(source))]
    else:
        labels = [tuple(label) for label in labels]
    if len(labels) != len(actual) or any(len(label)!=3 or any(type(v) is not int for v in label) for label in labels) or set(labels) != expected_labels:
        raise ValueError('labels must cover each copy and base index exactly once')
    max_distance = F(0)
    for (x,y,w),(a,b,i) in zip(actual,labels):
        sx,sy,sw = source[i]
        if w != sw:
            raise ValueError('matched mass differs from base')
        dx=(x-(sx+a*bx))%px; dy=(y-(sy+b*by))%py
        distance=max(min(dx,px-dx),min(dy,py-dy))
        if distance > delta:
            raise ValueError('matched displacement exceeds delta')
        max_distance=max(max_distance,distance)
    if sum((w for _,_,w in actual),F(0)) != nxcopy*nycopy*total:
        raise ValueError('total mass mismatch')
    bound=total/(bx*by)*(1+(bx+2*delta)/(2*tail))*(1+(by+2*delta)/(2*tail))
    if bound >= lam:
        raise ValueError('transport tail bound does not prove the requested radius cutoff')
    def canonical_digest(obj):
        return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    matching={'base_input':base,'base_canonical_sha256':canonical_digest(base),'data_canonical_sha256':canonical_digest(data),'copy_factors':[nxcopy,nycopy],'labels':labels,'rows':len(actual),'delta':str(delta),'maximum_torus_distance':str(max_distance),'base_mass':str(total),'expected_mass':str(nxcopy*nycopy*total),'base_periods':[str(bx),str(by)],'level':str(lam),'tail_radius':str(tail),'transport_tail_bound':str(bound),'original_code_sha256':'b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'}
    px, py, lam = (F(data[k]) for k in ('Px','Py','lam'))
    if px <= 0 or py <= 0 or lam <= 0 or tail <= 0 or bits < 0:
        raise ValueError('positive periods, level, radius and nonnegative bits required')
    raw_pts = [(F(x), F(y), F(w)) for x,y,w in data['pts']]
    if any(w < 0 for _,_,w in raw_pts):
        raise ValueError('negative weights are not allowed')
    pts = [(x%px, y%py, w) for x,y,w in raw_pts if w > 0]
    if not pts:
        raise ValueError('positive total mass required')
    mass = sum((w for _,_,w in pts), F(0))
    assert mass == nxcopy*nycopy*total
    D = 1 << bits
    # A shared integer scale makes every rectangle coordinate and area exact.
    scale = math.lcm(D, px.denominator, py.denominator, tail.denominator, *(v.denominator for p in pts for v in p[:2]))
    W = math.lcm(*(w.denominator for _,_,w in pts))
    R, X, Y = int(tail*scale), int(px*scale), int(py*scale)
    atoms = {}
    for x,y,w in pts:
        xi, yi, wi = int(x*scale), int(y*scale), int(w*W)
        for i in range((-R-xi)//X, (X+R-xi)//X+1):
            xx=xi+i*X
            if not -R <= xx <= X+R: continue
            for j in range((-R-yi)//Y, (Y+R-yi)//Y+1):
                yy=yi+j*Y
                if -R <= yy <= Y+R: atoms[xx,yy]=atoms.get((xx,yy),0)+wi
    xs=sorted({x for x,y in atoms}); ys=sorted({y for x,y in atoms})
    nx,ny=len(xs),len(ys)
    p=[[0]*(ny+1) for _ in range(nx+1)]
    xm={x:i+1 for i,x in enumerate(xs)}; ym={y:i+1 for i,y in enumerate(ys)}
    for (x,y),w in atoms.items(): p[xm[x]][ym[y]] += w
    for i in range(1,nx+1):
        row,prev=p[i],p[i-1]; acc=0
        for j in range(1,ny+1):
            acc+=row[j]; row[j]=prev[j]+acc
    inner=set(); outer=set()
    count={'boxes':0,'zero_mass':0,'span_rejected':0,'root_checked':0,'inner_empty':0,'outer_empty':0}
    # s^2 = M/(4 W lambda); compare exact integer squares before sqrt.
    A=lam.denominator; B=4*W*lam.numerator
    for i,a in enumerate(xs):
      for k in range(i,bisect.bisect_right(xs,a+2*R)):
        b=xs[k]
        if b-R >= X or a+R <= 0: continue
        for j,c in enumerate(ys):
          for l in range(j,bisect.bisect_right(ys,c+2*R)):
            d=ys[l]
            if d-R >= Y or c+R <= 0: continue
            count['boxes']+=1
            M=p[k+1][l+1]-p[i][l+1]-p[k+1][j]+p[i][j]
            if not M: count['zero_mass']+=1; continue
            span=max(b-a,d-c)
            if span*span*B >= 4*scale*scale*M*A:
                count['span_rejected']+=1; continue
            q=M*A*D*D//B
            root=math.isqrt(q)
            assert root*root*B <= M*A*D*D < (root+1)*(root+1)*B
            count['root_checked']+=1
            lo=min(R,max(0,root-1)*(scale//D)); hi=min(R,(root+1)*(scale//D))
            for s, target, label in ((lo,inner,'inner_empty'),(hi,outer,'outer_empty')):
                x0,x1=max(0,b-s),min(X,a+s)
                y0,y1=max(0,d-s),min(Y,c+s)
                if x0<x1 and y0<y1: target.add((x0,x1,y0,y1))
                else: count[label]+=1
    ai,ao=union_area(inner),union_area(outer)
    denom=scale*scale
    result={'period_x':str(px),'period_y':str(py),'level':str(lam),'mass':str(mass),'tail_radius':str(tail),'tail_bound':str(bound),'tail_bound_float':float(bound),'bits':bits,'scale':str(scale),'weight_scale':str(W),'input_atoms':len(pts),'window_atoms_distinct':len(atoms),'window_x_coordinates':nx,'window_y_coordinates':ny,'counts':count,'inner_rectangles_unique':len(inner),'outer_rectangles_unique':len(outer),'inner_area':str(F(ai,denom)),'outer_area':str(F(ao,denom)),'inner_ratio':str(F(ai,denom)*lam/mass),'outer_ratio':str(F(ao,denom)*lam/mass),'inner_ratio_float':float(F(ai,denom)*lam/mass),'outer_ratio_float':float(F(ao,denom)*lam/mass),'python':sys.version,'platform':platform.platform()}
    result['transport_matching']=matching
    if rectangles_out:
        with open(rectangles_out,'w') as f: json.dump({'scale':str(scale),'inner':[list(map(str,r)) for r in sorted(inner)],'outer':[list(map(str,r)) for r in sorted(outer)]},f,separators=(',',':'))
    return result


def main():
    p=argparse.ArgumentParser(); p.add_argument('input'); p.add_argument('--base',required=True); p.add_argument('--copy-factors',nargs=2,type=int,required=True); p.add_argument('--delta',required=True); p.add_argument('--labels'); p.add_argument('--bits',type=int,default=40); p.add_argument('--tail-radius',default='5'); p.add_argument('--output',required=True); p.add_argument('--rectangles-out'); args=p.parse_args()
    start=time.time(); data=json.loads(Path(args.input).read_text()); base=json.loads(Path(args.base).read_text()); labels=json.loads(Path(args.labels).read_text()) if args.labels else None
    result=calculate_transport(data,base,args.copy_factors,F(args.delta),args.bits,F(args.tail_radius),labels,args.rectangles_out)
    result.update(input_sha256=sha(args.input),base_sha256=sha(args.base),labels_sha256=sha(args.labels) if args.labels else None,code_sha256=sha(__file__),elapsed_seconds=time.time()-start)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
if __name__=='__main__': main()
