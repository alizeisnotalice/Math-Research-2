"""Inner polygon for arbitrary independent translations and rotations in one eleven-coordinate box."""
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
from hashlib import sha256
import sys
import numpy as np
from scipy.spatial import HalfspaceIntersection
from hinge_bank import S,edges,det,q,write_poly
from flint import arb,ctx

def core(comp,box,k=128):
    assert len(box)==11 and all(type(a) is int and type(b) is int and a<=b for a,b in box)
    assert sha256((comp/'certify_nine.py').read_bytes()).hexdigest()=='610e94cf22e901453bc34f77afb8f0994d9e471432cbe20498fca3b1530c9541'
    sys.path.insert(0,str(comp.resolve()))
    from certify_nine import hull,append_point,arc_body,S as old
    ctx.prec=160;G=2**44;pi=arb.pi()
    disk=[]
    for i in range(8*k):
        a=pi*i/(4*k);append_point(disk,a.cos()/2,a.sin()/2)
    D=old+4
    # Produce a small integer lower bound in each sampled support direction.
    local={}
    for n,j in ((3,None),(5,2),(7,5),(9,8)):
        for idx in range(1 if j is None else 2):
            v=[];ang=0 if j is None else F(box[j][idx],G)
            arc_body(v,k,pi,n,q(ang),arb(0),arb(0));local[n,idx]=v
    # A centered regular Reuleaux polygon has continuous support derivative;
    # h'' = -h on a vertex sector, and h''=1-h on a circular-arc sector.
    # In particular |h''|<=1. Every intermediate support exceeds the lesser
    # endpoint support minus width^2/8. Normalize support on the unit circle.
    eps=max(F((box[j][1]-box[j][0])**2,8*G*G) for j in (2,5,8))
    # Including independent rounding of a moved center, all support errors are
    # at most 3/old. Working normal inequalities on raw /old coordinates.
    err=eps+F(3,old)
    constraints={}
    def emit(a,b,c):
        # x is in the final 1/S grid. Directional corner hull has support
        # c/old when a,b is an edge normal in its original integer coordinates.
        m=isqrt(a*a+b*b)+1
        h=(F(S*c,old)-S*m*err).__floor__()
        key=(a,b)
        if key not in constraints or constraints[key]>h:constraints[key]=h
    for sx,sy in product((1,-1),repeat=2):
        choices={}
        for n,j in ((3,None),(5,2),(7,5),(9,8)):
            t=0 if j is None else j+1
            x=round(F(box[t][int(sx<0)]*old,G));y=round(F(box[t+1][int(sy<0)]*old,G))
            for z in range(1 if j is None else 2):
                choices[n,z]=[(a+x,b+y) for a,b in local[n,z]]
        base=disk+choices[3,0]
        for a,b,c in product((0,1),repeat=3):
            h=hull(base+choices[5,a]+choices[7,b]+choices[9,c])
            for p,r in edges(h):
                x,y=r[1]-p[1],p[0]-r[0]
                if x*sx>=0 and y*sy>=0:emit(x,y,x*p[0]+y*p[1])
            # Cardinal rays bound the angular ends even when no edge is axial.
            for x,y in ((sx,0),(0,sy)):
                emit(x,y,max(x*p+y*v for p,v in h))
    lines=[(a,b,c) for (a,b),c in constraints.items()]
    mat=np.array([(a,b,-c/S) for a,b,c in lines],float)
    mat/=np.linalg.norm(mat[:,:2],axis=1)[:,None]
    verts=HalfspaceIntersection(mat,np.zeros(2)).intersections
    # Inward rounding of the numerically suggested polygon; every resulting
    # integer vertex is checked against every exact necessary half-plane.
    poly=hull([(round(float(x)*(S-7)),round(float(y)*(S-7))) for x,y in verts])
    assert all(a*x+b*y<=c for x,y in poly for a,b,c in lines)
    return poly

if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    for x in ('comp','boxes','output'):p.add_argument(x,type=Path)
    p.add_argument('--k',default=128,type=int);a=p.parse_args();assert __debug__
    b=json.loads(a.boxes.read_bytes())
    with a.output.open('w') as f:
        f.write(f'{S} {len(b)}\n')
        for row in b:write_poly(f,core(a.comp,row,a.k))
