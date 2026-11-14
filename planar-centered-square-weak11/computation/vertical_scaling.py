#!/usr/bin/env python3
"""Exact rectangle enclosures for a vertical scaling family of periodic halos.

All source boxes come from one enlarged periodic window.  Their dyadic inner
and outer radii do not depend on beta.  interval_upper returns a superset of
every outer halo in a parameter interval, including the moving cell boundary.
The cutoff is valid only when the checked uniform counting bound is below 1.
"""
import argparse
import bisect
import hashlib
import json
import math
import platform
import sys
import time
from fractions import Fraction as F
from pathlib import Path

from periodic_halo import union_area


class VerticalFamily:
    def __init__(self, data, beta_min=F(999, 1000), beta_max=F(1001, 1000), bits=40, radius=F(5)):
        self.bmin, self.bmax, self.radius = F(beta_min), F(beta_max), F(radius)
        self.px, self.py, self.lam = (F(data[k]) for k in ('Px', 'Py', 'lam'))
        if not (0 < self.bmin <= self.bmax and min(self.px, self.py, self.lam, self.radius) > 0 and bits >= 0):
            raise ValueError('invalid family parameters')
        pts = [(F(x) % self.px, F(y) % self.py, F(w)) for x, y, w in data['pts']]
        if any(w < 0 for _, _, w in pts):
            raise ValueError('negative mass')
        pts = [(x, y, w) for x, y, w in pts if w]
        self.mass = sum((w for _, _, w in pts), F(0))
        if not self.mass:
            raise ValueError('zero mass')
        # g(R,beta) decreases in beta: this establishes the cutoff uniformly.
        self.tail_bound = self.mass / (self.px*self.bmin*self.py) * (1+self.px/(2*self.radius)) * (1+self.bmin*self.py/(2*self.radius))
        if self.tail_bound >= self.lam:
            raise ValueError('uniform tail bound fails')
        D = 1 << bits
        S = math.lcm(D, self.px.denominator, self.py.denominator, self.radius.denominator, *(v.denominator for p in pts for v in p[:2]))
        W = math.lcm(*(w.denominator for _, _, w in pts))
        self.scale, self.bits = S, bits
        X, Y, R = int(self.px*S), int(self.py*S), int(self.radius*S)
        atoms = {}
        # An expanded nominal y window covers the actual window for all beta.
        margin = self.radius/self.bmin*S
        for x, y, w in pts:
            xi, yi, wi = int(x*S), int(y*S), int(w*W)
            for i in range((-R-xi)//X, (X+R-xi)//X+1):
                xx = xi+i*X
                if not -R <= xx <= X+R:
                    continue
                for j in range(math.floor((-margin-yi)/Y), math.floor((Y+margin-yi)/Y)+1):
                    yy = yi+j*Y
                    if -margin <= yy <= Y+margin:
                        atoms[xx, yy] = atoms.get((xx, yy), 0)+wi
        xs, ys = sorted({x for x,y in atoms}), sorted({y for x,y in atoms})
        p = [[0]*(len(ys)+1) for _ in range(len(xs)+1)]
        xm, ym = {x:i+1 for i,x in enumerate(xs)}, {y:i+1 for i,y in enumerate(ys)}
        for (x,y), w in atoms.items():
            p[xm[x]][ym[y]] += w
        for i in range(1,len(xs)+1):
            acc = 0
            for j in range(1,len(ys)+1):
                acc += p[i][j]
                p[i][j] = p[i-1][j]+acc
        self.boxes = []
        count = 0
        A, B = self.lam.denominator, 4*W*self.lam.numerator
        # Only necessary conditions that hold across the entire family are used
        # to remove boxes, always with the outer radius or the cutoff R.
        for i, a in enumerate(xs):
            for k in range(i, bisect.bisect_right(xs, a+2*R)):
                b = xs[k]
                if b-R >= X or a+R <= 0:
                    continue
                for j, c in enumerate(ys):
                    for l in range(j, bisect.bisect_right(ys, c+int(2*R/self.bmin))):
                        d = ys[l]
                        if self.bmin*(d-c) >= 2*R:
                            continue
                        if min(self.bmin*d, self.bmax*d)-R >= self.bmax*Y or max(self.bmin*c, self.bmax*c)+R <= 0:
                            continue
                        count += 1
                        M = p[k+1][l+1]-p[i][l+1]-p[k+1][j]+p[i][j]
                        if not M:
                            continue
                        q = M*A*D*D//B
                        root = math.isqrt(q)
                        assert root*root*B <= M*A*D*D < (root+1)*(root+1)*B
                        lo = min(R, max(0, root-1)*(S//D))
                        hi = min(R, (root+1)*(S//D))
                        if b-a >= 2*hi or self.bmin*(d-c) >= 2*hi:
                            continue
                        # Do not discard boxes whose y receiver is currently empty.
                        self.boxes.append((a,b,c,d,lo,hi))
        self.enumeration = {'window_atoms_distinct':len(atoms), 'x_coordinates':len(xs), 'y_coordinates':len(ys), 'candidate_boxes':count, 'retained_boxes':len(self.boxes)}

    def _rectangles(self, lo, hi, inner=False):
        lo, hi = F(lo), F(hi)
        if not self.bmin <= lo <= hi <= self.bmax:
            raise ValueError('beta interval outside the certified family')
        if inner and lo != hi:
            raise ValueError('inner rectangles require a single beta')
        Q = math.lcm(lo.denominator, hi.denominator)
        L, H = int(lo*Q), int(hi*Q)
        X, Y = int(self.px*self.scale), int(self.py*self.scale)
        rects = set()
        for a,b,c,d,rlo,rhi in self.boxes:
            r = rlo if inner else rhi
            x0, x1 = max(0,b-r)*Q, min(X,a+r)*Q
            # Each moving endpoint is affine, hence its minimum/maximum occurs
            # at an interval endpoint. Clipping to max cell gives a superset.
            y0 = max(0, min(L*d,H*d)-r*Q)
            y1 = min(H*Y, max(L*c,H*c)+r*Q)
            if x0 < x1 and y0 < y1:
                rects.add((x0,x1,y0,y1))
        return rects, self.scale*Q

    def evaluate(self, beta):
        beta = F(beta)
        inner, si = self._rectangles(beta,beta,True)
        outer, so = self._rectangles(beta,beta)
        assert si == so
        return F(union_area(inner),si*si), F(union_area(outer),so*so)

    def interval_upper(self, lo, hi):
        rects, s = self._rectangles(lo,hi)
        return F(union_area(rects),s*s)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('--betas', nargs='*', default=['1'])
    parser.add_argument('--interval', nargs=2)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = time.time()
    input_path = Path(args.input)
    family = VerticalFamily(json.loads(input_path.read_text()))
    out = {'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(), 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'python':sys.version, 'platform':platform.platform(), 'bits':family.bits, 'beta_min':str(family.bmin), 'beta_max':str(family.bmax), 'tail_bound':str(family.tail_bound), 'enumeration':family.enumeration, 'evaluations':{}}
    for arg in args.betas:
        lo,hi = family.evaluate(F(arg))
        out['evaluations'][arg] = {'inner_area':str(lo), 'outer_area':str(hi), 'inner_ratio':str(lo*family.lam/family.mass), 'outer_ratio':str(hi*family.lam/family.mass), 'inner_ratio_float':float(lo*family.lam/family.mass), 'outer_ratio_float':float(hi*family.lam/family.mass)}
    if args.interval:
        lo,hi = map(F,args.interval)
        upper = family.interval_upper(lo,hi)
        out['interval'] = {'lo':str(lo), 'hi':str(hi), 'upper_area':str(upper), 'upper_ratio':str(upper*family.lam/family.mass), 'upper_ratio_float':float(upper*family.lam/family.mass)}
    out['elapsed_seconds'] = time.time()-started
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
