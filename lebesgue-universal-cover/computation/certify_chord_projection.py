"""Uniform cap-area/parallel-chord support cuts for the complete moving 8D base box."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
import json
from math import ceil,floor,gcd,isqrt,sqrt
from pathlib import Path

from flint import arb,ctx
import numpy as np
from certify_illumination_projection import CENTER,GRID,TAU,inner_polygon,certify_rows
from certify_three_body_refinement import CHECKER_SHA

HALF=Q(1,4096)

def sharp_area_lower(source):
    path=Path(source)/'checker/lbcheck.py';assert sha256(path.read_bytes()).hexdigest()==CHECKER_SHA
    spec=spec_from_file_location('tomar',path);lib=module_from_spec(spec);spec.loader.exec_module(lib);ctx.prec=112
    ck=lib.Checker([3,5,7],float(TAU),768,3e-7);c=np.array(list(map(float,CENTER)));w=HALF+Q(1,2**50);h=np.full(8,float(w))
    for x,xc in zip(CENTER,c):assert abs(x-Q(float(xc)))+HALF<=Q(float(w))==w
    runs,pa=ck.verified_runs(c,h,c);assert runs is not None
    f,g,e=ck.taylor_model(runs,pa,h);l=f-e-sum(abs(z)*arb(float(w)) for z in g);n=floor(float(l.lower())*10**12)-1;assert l>arb(n)/10**12
    return Q(n,10**12)

def slices(p,h,b,off,levels):
    candidates=[]
    for dd in levels:
        level=max(h)-dd*off;borders=[]
        for i,vv in enumerate(h):
            j=(i+1)%len(h)
            if (vv>=level)==(h[j]>=level):continue
            bary=(level-h[j])/(vv-h[j]);borders.append(bary*p[i]+(1-bary)*p[j])
        if len(borders)!=2:continue
        ell=abs(borders[0]-borders[1]);assert ell>0 and level<=b
        candidates.append((level,ell))
    return candidates

def run(source,M):
    lower=sharp_area_lower(source);p,alpha=inner_polygon(HALF);points=[(Q(x,GRID),Q(y,GRID)) for x,y in p]
    delta=(Q(99,70)+Q(3,5))*HALF;outer=delta+Q(2,GRID)+Q(22,7)**2/(72*64**2)
    levels=list(map(Q,['.00001','.0001','.0005','.001','.002','.004','.008','.016','.032','.064']))
    cuts=[]
    for a,b in zip(p,p[1:]+p[:1]):
        nx=b[1]-a[1];ny=a[0]-b[0];g=gcd(abs(nx),abs(ny));nx//=g;ny//=g
        square=nx*nx+ny*ny;off=isqrt(square);off+=off*off<square;assert off**2>=square
        raw=[nx*x+ny*y for x,y in points];b=max(raw)+outer*off
        h=[alpha*z for z in raw];t=[alpha*(-ny*x+nx*y) for x,y in points];d=(TAU-lower)*square
        best=None
        for z,ell in slices(t,h,b,off,levels):
            dd=float(d);ee=float(ell);s=(dd+sqrt(dd*dd+2*dd*ee*float(b-z)))/ee
            # A float only proposes a rational cap intercept; its sufficiency is independently checked.
            bump=Q(ceil((s/off+1e-12)*10**12),10**12)*off
            upper=b+bump;assert bump>0 and ell*bump**2>2*d*(upper-z)
            if best is None or upper<best:best=upper
        assert best is not None;cuts.append(((nx,ny),best))
    rows,least,material=certify_rows(cuts,M)
    return {'tau':str(TAU),'center':[str(q) for q in CENTER],'box_halfwidth':str(HALF),'base_hull_area_lower':str(lower),'alpha':str(alpha),
            'angular_cells':M,'support_planes':len(cuts),'B9_material_points':len(material),'minimum_dual_gap_lower':str(Q(least,10**12)),
            'rows':rows,'scope':'strict local 8D box for every B9 orientation and translation; uses certified moving-base area independently of common inner polygon'}

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('source');ap.add_argument('--cells',type=int,default=256);a=ap.parse_args();print(json.dumps(run(a.source,a.cells),separators=(',',':')))
