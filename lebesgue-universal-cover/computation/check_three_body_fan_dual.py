"""Two-fan rational dual over all translation corners and an angle interval."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from itertools import combinations
import json
from math import floor
from pathlib import Path

from flint import arb
import numpy as np
from certify_three_body_refinement import CHECKER_SHA,corner_bounds

SCALE=10**12

def pair_max(a,b):
    slopes=[x-y for x,y in zip(a,b)]
    candidates={Q(0),Q(1)}
    for i,j in combinations(range(len(a)),2):
        if slopes[i]!=slopes[j]:
            v=Q(b[j]-b[i], slopes[i]-slopes[j]);
            if 0<v<1: candidates.add(v)
    return max((min(v*m+n for m,n in zip(slopes,b)),v) for v in candidates)


def demo(source):
    path=Path(source)/'checker/lbcheck.py';assert sha256(path.read_bytes()).hexdigest()==CHECKER_SHA
    spec=spec_from_file_location('tomar',path);lib=module_from_spec(spec);spec.loader.exec_module(lib)
    # Exact dyadic argument: (t3x,t3y,angle5,t5x,t5y).
    center=[-10,-8,116,102,-13];half=[6,6,18,6,6];grid=4096
    c=np.array(center)/grid; h=np.array(half)/grid
    checker=lib.Checker([3,5],0.83475,256,1e-5)
    rows={}
    shifts=[(0,0,0,0,0)]+[tuple(s if i==j else 0 for i in range(5)) for j in range(5) for s in (-1,1)]
    for shift in shifts:
        runs,pa=checker.verified_runs(c,h,c+h*np.array(shift))
        if runs is None: continue
        corners=[]
        for ball in corner_bounds(runs,pa,h):
            v=floor(float(ball.lower())*SCALE)-1
            assert ball>arb(v)/SCALE
            corners.append(v)
        rows[shift]=corners
    single=max(min(row) for row in rows.values())
    pairs=[]
    for (sa,a),(sb,b) in combinations(rows.items(),2):
        val,wt=pair_max(a,b);pairs.append((val,wt,sa,sb))
    (dual,w,sa,sb)=max(pairs)
    assert 0<w<1 and dual>single
    a,b=rows[sa],rows[sb]; bound=min((w*x+(1-w)*y)/SCALE for x,y in zip(a,b)); assert bound==dual/SCALE
    return {'grid':grid,'center':center,'halfwidth':half,'shifts':[sa,sb],
            'weights':[str(w),str(1-w)],'corner_lower_integer_rows':[a,b],'scale':SCALE,
            'best_single_among_tested':str(Q(single,SCALE)),'combined_lower':str(bound),
            'combined_display':float(bound),'gain_display':float((dual-single)/SCALE),'valid_fans':len(rows),
            'scope':'one closed five-dimensional dyadic box; all angle and translation coordinates included'}

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('source');args=ap.parse_args();print(json.dumps(demo(args.source),indent=2))
