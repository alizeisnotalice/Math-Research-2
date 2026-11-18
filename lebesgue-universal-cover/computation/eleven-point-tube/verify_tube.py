"""Replay the full line segment of independent pose parameters and rigorously thicken it."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import gzip
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
from template import build
from verify_angles import rational, verify
from flint import arb,ctx


def check_one(job):
    j,N,A,k,K,row,comp,points=job
    t=F(2*j+1,2*N)
    if row==0:
        h,d=build(comp,False,k,True,fraction=t)
        twice=sum(h[i-1][0]*h[i][1]-h[i-1][1]*h[i][0] for i in range(len(h)))
        floor=F(twice,2*d*d)
        assert floor>A
        return (0,0,floor)
    assert isinstance(row,list) and row
    data=dict(seven=False,five_only=True,target=str(A),samples=k,fraction=str(t),period=K,cells=row)
    result=verify(comp,points,data)
    return (1,result['cases'],None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('computation','points','certificate','out'):parser.add_argument(key,type=Path)
    parser.add_argument('--workers',default=2,type=int)
    a=parser.parse_args()
    assert __debug__ and a.workers>=1
    ctx.prec=160
    data=json.loads(gzip.decompress(a.certificate.read_bytes()))
    N=data['N'];A=F(data['target']);k=data['samples'];K=data['period']
    assert (N,k,A,K)==(1024,128,F('.837595'),131072) and len(data['cells'])==N
    inputs=['free-p30-reoptimized-pose.json','seven-point-families-pose.json','certify_nine.py']
    src=[(a.computation/z).read_bytes() for z in inputs]
    base=list(map(F,json.loads(src[0])))
    end=list(map(F,json.loads(src[1])['pose']))
    L={3:F('.289'),5:F('.163'),7:F('.115'),9:F('.089')}
    # Tighter constants need the Arb bounds themselves for the movement along the segment.
    speeds=[]
    for n,j in ((3,None),(5,2),(7,5),(9,8)):
        pos=0 if j is None else j+1
        speed=rational((end[pos]-base[pos])**2+(end[pos+1]-base[pos+1])**2).sqrt()
        angular=(arb.pi()/(2*n)).tan()/2
        assert angular<rational(L[n])
        if j is not None:speed+=angular*rational(abs(end[j]-base[j]))
        assert speed<rational(F('.0201'))
        speeds.append(str(speed))
    with ProcessPoolExecutor(a.workers) as pool:
        rs=list(pool.map(check_one,((i,N,A,k,K,r,a.computation,a.points) for i,r in enumerate(data['cells'])),chunksize=4))
    chords=sum(r[0] for r in rs);direct=[r[2] for r in rs if not r[0]]
    assert len(direct)+chords==N and 0<chords<N
    scatter=F('.0201')/(2*N)
    margins={}
    for target,extra in ((F('.83755'),F('.00000361')),(F('.83750'),F('.00001854'))):
        lower=A/(1+2*(scatter+extra))**2
        assert lower>target
        margins[str(target)]=dict(independent_motion_extra=str(extra),area_strict_lower=str(lower))
    out=dict(cells=N,chord_cells=chords,curve_only_cells=len(direct),
             angular_leaf_checks=sum(r[1] for r in rs),minimum_curve_only_midpoint=str(min(direct)),
             pose_speed_enclosures=speeds,subinterval_motion_budget=str(scatter),
             midpoint_strict_lower=str(A),uniform_tube=margins,arc_samples=k,python=platform.python_version(),
             certificate_sha256=sha256(a.certificate.read_bytes()).hexdigest(),
             point_sha256=sha256(a.points.read_bytes()).hexdigest(),
             source_sha256={n:sha256(b).hexdigest() for n,b in zip(inputs,src)})
    a.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
