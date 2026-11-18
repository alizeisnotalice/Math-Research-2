"""Propose Q11 angle rows on sampled genuine Tomar source boxes."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from certify_finite_point_leaf import U
from certify_tomar_six_high import SourceTree

EPS = Fraction(1,10**9)

def reflected_average(d):
    cloud = d["points"]
    assert len(cloud)==11 and cloud[0]==[0,0]
    den = d["point_denominator"]
    assert type(den) is int and den>0
    assert all(len(z)==2 and all(type(t) is int for t in z) for z in cloud)
    assert max((x-X)**2+(y-Y)**2 for (x,y),(X,Y) in combinations(cloud,2))<den**2
    pairs = [(x+X,y-Y) for (x,y),(X,Y) in zip(cloud,cloud[:1]+cloud[:0:-1],strict=True)]
    assert all(Fraction((a-2*x)**2+(b-2*y)**2,4*den**2)<EPS**2
               for (a,b),(x,y) in zip(pairs,cloud,strict=True))
    return [list(z) for z in pairs],2*den

def main():
    from generate_multi_host_points import Sampler, m, u
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("tree",type=Path)
    p.add_argument("q11",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("cases",nargs="+")
    a=p.parse_args()
    tree=SourceTree(a.tree)
    keys=[tuple(map(int,key.split("/"))) for key in a.cases]
    centers=np.array([[float((x+y)/2) for x,y in tree.leaf(*key)] for key in keys])
    cloud,den=reflected_average(json.loads(a.q11.read_bytes()))
    targets=np.asarray(cloud,dtype=float)/den
    normals=[]
    normal_ids={}
    def index(z):
        d=tuple(map(int,np.rint(z*U)))
        if d not in normal_ids:
            normal_ids[d]=len(normals)
            normals.append(d)
        return normal_ids[d]
    sampler=Sampler(centers,384,360)
    rows=[[] for _ in range(384)]
    for H,sign,i,theta,jj,mat,weight,gap,xy in sampler.angles_rows(targets):
        if sign!=1:
            continue
        normal=sampler.host_table[H][0]
        v=(*map(index,normal[jj]),*map(int,mat))
        if v not in rows[i]:
            rows[i].append(v)
    tris=u(2*np.pi*np.arange(3)/3)
    fallback=list(map(index,tris))
    for i,theta in enumerate(sampler.angles):
        if not rows[i]:
            rows[i].append((*fallback,*map(int,np.argmax((tris@m(theta))@targets.T,axis=1))))
    a.out.write_text(json.dumps({"cases":keys,"directions":normals,"cells":rows},separators=(",",":"))+"\n")
    print(len(normals),sum(map(len,rows)),flush=True)

if __name__=="__main__":
    main()
