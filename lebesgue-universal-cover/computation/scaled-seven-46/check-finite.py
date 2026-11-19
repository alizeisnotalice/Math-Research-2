"""Verify every side of five exported rational polygons and their diameters."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip
import json

from flint import arb, ctx

POINTS_SHA = "5e39cccdc5ee67c8e32cb3d1112c1a2caac89717289cdc24d84411fc3ec82c3e"
SHAPES = (0, 3, 5, 7, 9)

def cross(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])

def polygon(points):
    points=sorted(set(points))
    def chain(rows):
        h=[]
        for a in rows:
            while len(h)>1 and cross(h[-2],h[-1],a)<=0:
                h.pop()
            h.append(a)
        return h[:-1]
    return chain(points)+chain(reversed(points))

def distance2(a,b):
    return (a[0]-b[0])**2+(a[1]-b[1])**2

def rotating_diameter(h):
    size=len(h); j=1; best=0
    assert size>=3
    for i in range(size):
        a,b=h[i],h[(i+1)%size]
        while cross(a,b,h[(j+1)%size])>cross(a,b,h[j%size]):
            j+=1
        opp=h[j%size]
        best=max(best,distance2(a,opp),distance2(b,opp))
        nxt=h[(j+1)%size]
        if cross(a,b,nxt)==cross(a,b,opp):
            best=max(best,distance2(a,nxt),distance2(b,nxt))
    assert j<2*size+1
    return best

def support_geometry(n):
    pi=arb.pi()
    axes=[]
    for k in range(n):
        s,c=(2*k*pi/n).sin_cos()
        axes.append((c,s))
    half=(pi/(2*n)).sin_cos()
    return axes,half[0],1/(2*half[1]),*(pi/n).sin_cos()

def height_upper(nx,ny,length,geom):
    if geom is None:
        return length/2
    axes,transition,radius,sin_step,cos_step=geom
    projs=[nx*c+ny*s for c,s in axes]
    best=max(x.lower() for x in projs)
    upper=[]
    for (c,s),dot in zip(axes,projs,strict=True):
        if dot.upper()<best:
            continue
        transverse=abs(ny*c-nx*s)
        boundary=transition*length
        if not transverse>boundary:
            upper.append((radius*dot).upper())
        if not transverse<boundary:
            upper.append((length-radius*(cos_step*dot+sin_step*transverse)).upper())
    return max(upper)

def check(directory):
    datafile=directory/'rational-seven-83746.json.gz'
    payload=gzip.decompress(datafile.read_bytes())
    assert sha256(payload).hexdigest()==POINTS_SHA
    groups=json.loads(payload)
    report=json.loads((directory/'report.json').read_bytes())
    assert len(groups)==7 and report['normalized_area_lower']=='18415940155945779266357421875/21990232335617674994555813888'
    factors=[];summary=[]
    for cloud,rec,n in zip(groups[:5],report['rows'][:5],SHAPES,strict=True):
        a,b=rec['integer_axis']
        d=cloud['denominator'];assert a*a+b*b==d
        pts=[]
        for x,y in cloud['points']:
            xx,yy=x*a-y*b,x*b+y*a
            assert xx%d==yy%d==0
            pts.append((xx//d,yy//d))
        h=polygon(pts)
        assert rotating_diameter(h)==d
        geom=support_geometry(n) if n else None
        minimum=None
        for u,v in zip(h,h[1:]+h[:1],strict=True):
            nx,ny=v[1]-u[1],u[0]-v[0]
            norm=arb(nx*nx+ny*ny).sqrt()
            z=height_upper(nx,ny,norm,geom)
            excess=((nx*u[0]+ny*u[1]-(10**21-2)*z)/norm).lower()
            assert excess>0,(n,u,v,str(excess))
            if minimum is None or excess<minimum:
                minimum=excess
        factors.append(Fraction((10**21-2)**2,d))
        summary.append(dict(n=n,vertices=len(h),integer_diameter_squared=str(d),minimum=str(minimum)))
    for cloud in groups[5:]:
        pts=[tuple(p) for p in cloud['points']]
        sq=max(distance2(a,b) for a,b in combinations(pts,2))
        assert sq==rotating_diameter(polygon(pts))<cloud['denominator']**2
    bound=min(factors)*Fraction(report['normalized_area_lower'])
    assert bound==Fraction(report['true_area_lower'])>Fraction('.837460000000005871546308')
    return dict(walls=sum(x['vertices'] for x in summary),shapes=summary,
                area_lower=str(bound),uncompressed_point_sha256=POINTS_SHA,precision=ctx.prec)

if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument('data',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();ctx.prec=192
    result=check(a.data)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['walls'],result['area_lower'])
