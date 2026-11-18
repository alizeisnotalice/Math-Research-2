"""Build rational polygons common to each signed low-area class of an original leaf."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import sys

from flint import arb,ctx

FINE=2**44
COORD=2**28
RADII=(F('.42'),F('.47'),F('.48'),F('.49'))
SLOPE=(F(0), F('.163'),F('.115'),F('.089'))


def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]


def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]


def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def polygon_area(poly):
    return sum(cross(a,b) for a,b in zip(poly,poly[1:]+poly[:1],strict=True))/2


def dual_polygon(poly):
    assert polygon_area(poly)>0
    planes=[]
    for p,q in zip(poly,poly[1:]+poly[:1],strict=True):
        normal=(q[1]-p[1],p[0]-q[0])
        height=dot(normal,p)
        assert height>0
        planes.append((normal[0]/height,normal[1]/height))
    return planes


def polar(body):
    out=[]
    for p,q in zip(body,body[1:]+body[:1],strict=True):
        z=cross(p,q)
        assert z>0
        out.append(((q[1]-p[1])/z,(p[0]-q[0])/z))
    return out


def integer_sqrt(n):
    s=isqrt(n)
    return s+(s*s<n)


def inner(bounds,high,table,steps,sizes,hull):
    lo,hi=bounds[:11],bounds[11:]
    points=[(F(int(x),COORD),F(int(y),COORD)) for x,y in table[0,:sizes[0]]]
    for b in range(4):
        j=3*b
        x,y=F(lo[j]+hi[j],2*FINE),F(lo[j+1]+hi[j+1],2*FINE)
        tw=F(integer_sqrt((hi[j]-lo[j])**2+(hi[j+1]-lo[j+1])**2),2*FINE)
        row=0
        if b:
            reference=high[j-1]
            row=((lo[j-1]+hi[j-1])*1024+reference)//(2*reference)
            assert 0<=row<=1024
            angle_error=F(max(abs(1024*lo[j-1]-row*reference),
                              abs(1024*hi[j-1]-row*reference)),1024*FINE)
            tw+=SLOPE[b]*angle_error
        scale=1-tw/RADII[b]
        assert scale>0
        for u,v in table[row,steps[b+1]:steps[b+2]]:
            points.append((x+scale*F(int(u),COORD),y+scale*F(int(v),COORD)))
    return hull(points)


def run(comp,tree,table,scan):
    ctx.prec=160
    sys.path.insert(0,str(comp))
    from certify_nine import hull
    from certify_tomar_six_high import SourceTree
    from verify_tomar_fan_table import verified_table
    data=read(scan)
    raw=SourceTree(tree)
    box=raw.leaf(*data['key'])
    for low,high in zip(*data['root'],strict=True):
        assert low<=high
    assert data['root']==[[int(x*FINE) for x in row] for row in zip(*box,strict=True)]
    table_values,sizes=verified_table(table,tree,raw.root)
    table_values=table_values.reshape((1025,sum(sizes),2))
    high=[int(F(t)*FINE) for t in raw.root['hi']]
    steps=[0]
    for m in sizes:
        steps.append(steps[-1]+int(m))
    pi=arb.pi()
    for n,r,l in zip((3,5,7,9),RADII,SLOPE,strict=True):
        assert 1-1/(2*(pi/(2*n)).cos())>ar(r)
        assert n==3 or (pi/(2*n)).tan()/2<ar(l)
    results=[]
    for group,slice_ in zip(data['live_boxes'],data['slices'],strict=True):
        assert len(group)==slice_[2]>0
        common_dual=[]
        for box in group:
            common_dual.extend(dual_polygon(inner(box,high,table_values,steps,sizes,hull)))
        p=polar(hull(common_dual))
        area=polygon_area(p)
        assert 0<area<F(data['target'])
        floor=int(area*10**12)/10**12
        results.append(dict(live_boxes=len(group),facets=len(p),area_floor=floor,
                            area_exact=str(area),polygon=[list(map(str,r)) for r in p]))
    return dict(classification_axis=data['classification_axis'],tree_key=data['key'],target=data['target'],
                source_scan_sha256=sha(scan),common_inner_certificates=results)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    for name in ('comp','tree','table','scan'):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.comp,a.tree,a.table,a.scan),separators=(',',':')))
