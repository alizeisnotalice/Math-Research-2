"""Uniform five-parameter box via centered rational material fans."""
from pathlib import Path
from fractions import Fraction as F
import contextlib, io, itertools, json, runpy
with contextlib.redirect_stdout(io.StringIO()):
    data=runpy.run_path(str(Path(__file__).with_name("certify_centered_three_body.py")))
disk,tri,pent,hull,den,r=(data[k] for k in ('disk','tri','pent','hull','den','r'))
cells=512
# For a fixed angular cell, rotating a vertex of radius < .53 incurs
# error <= .53 * .106 / cells. Add the two rounding stages.
err=F(53,100)*F(106,1000*cells)+F(4,den)
assert F(1,2)/data['trig'](data['scale'](data['PI'],F(1,10)),'cos')[0]<F(53,100)
def det(p,q):
    return p[0]*q[1]-p[1]*q[0]
def sweep(step):
    worst=None
    for i in range(cells):
        t=F(106*(2*i+1),2000*cells)
        c=(1-t*t)/(1+t*t);s=2*t/(1+t*t)
        pent_i=[(round(c*x-s*y),round(s*x+c*y)) for x,y in pent]
        labels={p:0 for p in disk}
        labels.update({p:1 for p in tri});labels.update({p:2 for p in pent_i})
        h=hull(disk+tri+pent_i)
        fixed=0;moving=[]
        for p,q in zip(h,h[1:]+h[:1]):
            a,b=labels[p],labels[q]
            if a or b: moving.append((p,q,a,b))
            else:
                area=det(p,q);assert area>0
                fixed+=area
        for x,y,u,v in itertools.product((-step,step),repeat=4):
            shifts=((0,0),(x,y),(u,v))
            total=fixed
            for p,q,a,b in moving:
                z=det((p[0]+shifts[a][0],p[1]+shifts[a][1]),
                      (q[0]+shifts[b][0],q[1]+shifts[b][1]))
                if z<=0:return None
                total+=z
            worst=total if worst is None else min(worst,total)
    return F(worst,2*den**2)*(1-err/r)**2
for denominator in (250,500,1000,2000,4000):
    assert den%denominator==0
    lb=sweep(den//denominator)
    if lb is not None and lb>F(419,500):
        print(json.dumps(dict(halfwidth=str(F(1,denominator)),cells=cells,
            lower=str(lb),display=float(lb),vertex_error=str(err)),indent=2))
        break
else:raise ValueError("No certified width")
