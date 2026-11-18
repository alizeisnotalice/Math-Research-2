"""A rational polygon common to every convex universal subset of the moving P host."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from certify_nine import S, append_point, hull
from certify_sigma_mandatory_arcs import anchors, degree, g, verify


def common_polygon(left,right,parts,rays):
    assert parts>0 and rays>=4
    source=verify(left,right)
    dirs=[g.unit(2*g.PI*j/rays) for j in range(rays)]
    radii=[g.arb(2) for _ in dirs]
    for j in range(parts):
        a,b=(left+(right-left)*i/parts for i in (j,j+1))
        g.SIGMA=degree(a).union(degree(b))
        groups=[]
        for label in ('TI','GK','QY'):
            row=[(p,1-g.dot(p,p)) for p in anchors(label)[0].values()]
            assert min(c for _,c in row)>0
            groups.append(row)
        for i,u in enumerate(dirs):
            fixed=g.HALF
            for row in groups:
                reaches=[]
                for p,c in row:
                    d=g.dot(u,p)
                    reaches.append((d+(d*d+c).sqrt()).lower())
                fixed=max(fixed,min(reaches))
            radii[i]=min(radii[i],fixed)
    grid=[]
    for u,r in zip(dirs,radii,strict=True):
        append_point(grid,*g.scale(r,u))
    poly=hull(grid)
    two=sum(x*v-y*u for (x,y),(u,v) in zip(poly,poly[1:]+poly[:1]))
    area=F(two,2*(S+4)**2)
    return {'sigma_degrees':list(map(str,(left,right))),'sigma_parts':parts,'rays':rays,
            'common_polygon_vertices':len(poly),'common_denominator':S+4,'integer_vertices':poly,
            'exact_area':str(area),'area_lower_9digits':str(F((area*10**9).__floor__(),10**9)),
            'full_root_rectangle_counts':{k:v['rectangle_count'] for k,v in source['cases'].items()},
            'arc_verifier_sha256':sha256(Path(__file__).with_name('certify_sigma_mandatory_arcs.py').read_bytes()).hexdigest()}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('low',type=F)
    p.add_argument('high',type=F)
    p.add_argument('parts',type=int)
    p.add_argument('rays',type=int)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(common_polygon(args.low,args.high,args.parts,args.rays),indent=2)+'\n')
