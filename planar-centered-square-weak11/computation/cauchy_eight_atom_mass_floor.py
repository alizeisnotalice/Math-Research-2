#!/usr/bin/env python3
"""Finite rational-rectangle dual set certifying a uniform mass floor for all theta."""
import os,sys,json,hashlib,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import flint
from flint import arb,ctx
ctx.prec=128
ROOT=Path(__file__).resolve().parents[1]
AXIS=[(Q(1),Q(0)),(Q(-1),Q(0)),(Q(0),Q(1)),(Q(0),Q(-1))]
DIAG=[(Q(a,2),Q(b,2)) for a in [-1,1] for b in [-1,1]]

def A(q):return arb(q.numerator)/q.denominator
@lru_cache(None)
def F(q):
    if not q:return arb(0)
    a=A(q);return a*a.log()-a
@lru_cache(None)
def H(x,y):
    if not x or not y:return arb(0)
    sign=(1 if x>0 else -1)*(1 if y>0 else -1)
    return sign*(F(abs(x)+abs(y))-F(abs(x))-F(abs(y)))

def integral(box,z):
    xl,xh,yl,yh=box;a,b=z
    return H(xh-a,yh-b)-H(xl-a,yh-b)-H(xh-a,yl-b)+H(xl-a,yl-b)

def pointW(x,y,atoms):
    return 1/(x+y)-sum((1/(abs(x-a)+abs(y-b)) for a,b in atoms),Q(0))/4

def candidate(N):
    ys=[Q(j,N) for j in range(N+1)]
    xs=[Q(j,N) for j in range(N+1)]
    for left in [1,2,4,8,16]:xs.extend(Q(left)+Q(left*j,N) for j in range(1,N+1))
    boxes=[];va=arb(0);vd=arb(0)
    for xl,xh in zip(xs,xs[1:]):
        mult=4 if xh<=1 else 8
        for yl,yh in zip(ys,ys[1:]):
            x=(xl+xh)/2;y=(yl+yh)/2
            wa=pointW(x,y,AXIS);wd=pointW(x,y,DIAG)
            if wa+3*wd<=0:continue # exact theta=1/4 midpoint selection; selection need not ensure pointwise sign
            box=(xl,xh,yl,yh)
            u=integral(box,(Q(0),Q(0)))
            ia=u-sum((integral(box,z) for z in AXIS),arb(0))/4
            id_=u-sum((integral(box,z) for z in DIAG),arb(0))/4
            va+=mult*ia;vd+=mult*id_
            boxes.append({'box':[str(v) for v in box],'multiplicity':mult})
    return va,vd,boxes

if __name__=='__main__':
    assert os.path.ismount('/Volumes/data')
    st=time.monotonic();trials=[]
    for N in [8,16,32]:
        va,vd,boxes=candidate(N);trials.append({'N':N,'rectangles':len(boxes),'axis_integral':str(va),'diagonal_integral':str(vd)})
        if va>A(Q(40,19)) and vd>A(Q(40,19)):break
    else:raise RuntimeError('No certified mass floor by bounded refinement')
    out={'status':'PASS_all_theta_mass_strictly_gt_40_over_19','python':sys.version,'flint':flint.__version__,'precision_bits':128,'workers':1,'theta_range':'all real theta in [0,1]','mass_threshold':'40/19','selection':'exact midpoint W_axis+3W_diag>0; finite test set only, no pointwise positivity assumption','regions':'central first-quadrant unit square reflected four ways; x>1,y in[0,1] arm reflected/swapped eight ways; x<=32','integral_formula':'signed primitive H(x,y)=sign(x)sign(y)*(F(abs(x)+abs(y))-F(abs(x))-F(abs(y))), F(t)=t log t-t','axis_integral_on_S':str(va),'diagonal_integral_on_S':str(vd),'proof':'For every theta, integral(W_theta)+ >= integral_S W_theta = theta*V_axis+(1-theta)*V_diag >40/19. Multiplicities expand disjoint rectangles up to null boundaries.','N':N,'rectangles':boxes,'trials':trials,'elapsed_seconds':time.monotonic()-st,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=ROOT/'computation/results/cauchy-eight-atom-mass-floor.json';dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='rectangles'},indent=2))
