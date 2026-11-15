#!/usr/bin/env python3
"""Single-core exact geometry checks; no exhaustive halo-area claim."""
import json, platform, hashlib
from fractions import Fraction as F
from pathlib import Path

def rect_area(a):
    return max(F(0),a[1]-a[0])*max(F(0),a[3]-a[2])
def meet(a,b):
    return (max(a[0],b[0]),min(a[1],b[1]),max(a[2],b[2]),min(a[3],b[3]))
def integral(pieces,q):
    return sum((rho*rect_area(meet(s,q)) for s,rho in pieces),F(0))
def core(b,t): return (b[0]+t,b[1]-t,b[2]+t,b[3]-t)
def incore(x,b,t):
    return b[0]+t<x[0]<b[1]-t and b[2]+t<x[1]<b[3]-t

def main():
    R=F(1)
    # Contains arbitrarily narrow children, a four-block intersection, and empty cores.
    xs=[F(0),F(1,10),F(2),F(5)]
    ys=[F(-3),F(0),F(1,7),F(4)]
    blocks=[(a,b,c,d) for a,b in zip(xs,xs[1:]) for c,d in zip(ys,ys[1:])]
    pieces=[(b,F(j+1,7)) for j,b in enumerate(blocks)]
    points=[(x,y) for x in [F(0),F(1,20),F(1,10),F(1),F(2),F(3),F(5)]
                      for y in [F(-3),F(-1),F(0),F(1,14),F(1,7),F(1),F(4)]]
    checks=0
    for x in points:
      for r in [F(1,20),F(1,7),F(1,2),F(1)]:
        q=(x[0]-r,x[0]+r,x[1]-r,x[1]+r)
        # Only near-seam recipients require collars. Other centers use the owner core.
        if any(incore(x,b,R) for b in blocks): continue
        full=integral(pieces,q)
        saved=F(0)
        for b,rho in pieces:
          inner=core(b,2*R)
          lost=rect_area(meet(meet(q,b),inner)) if rect_area(inner)>0 else F(0)
          assert lost==0, (x,r,b,lost)
          saved+=rho*(rect_area(meet(q,b))-lost)
        assert full==saved
        checks+=1
    # Two 2x2 densities have identical complete coordinate marginals.
    diag=[((F(0),F(1),F(0),F(1)),F(1)),((F(1),F(2),F(1),F(2)),F(1))]
    anti=[((F(0),F(1),F(1),F(2)),F(1)),((F(1),F(2),F(0),F(1)),F(1))]
    q=(F(0),F(1,2),F(0),F(1,2))
    assert integral(diag,q)==F(1,4) and integral(anti,q)==0
    # Sharpness: R=1, collar w=3/2, child [0,10]x[-10,10].
    # f and g differ by y-translation; same collar=0 and same strict-core halo area.
    # At x=(7/8,0), Q_1 contains all f and no g. Total mass 5 > 4.
    sf=(F(13,8),F(7,4),F(-1,16),F(1,16))
    sg=(sf[0],sf[1],F(63,16),F(65,16))
    density=F(5)/rect_area(sf)
    q=(F(-1,8),F(15,8),F(-1),F(1))
    af=integral([(sf,density)],q)/4
    ag=integral([(sg,density)],q)/4
    assert af==F(5,4) and ag==0
    # In the capped problem g has zero mass in EVERY Q_r(x), r<=1.
    assert sg[2]>1
    out={'python':platform.python_version(),'arithmetic':'fractions.Fraction',
         'random_seed':None,'tolerance':None,'single_process':True,
         'collar_exact_checks':checks,'blocks':len(blocks),
         'minimum_block_width':str(F(1,10)),
         'marginal_counterexample_integrals':[str(integral(diag,q)) for q in [(F(0),F(1,2),F(0),F(1,2))]]+[str(F(0))],
         'w_3over2_sharpness_averages':[str(af),str(ag)],
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).parent/'results/recursive_geometry_exact.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
