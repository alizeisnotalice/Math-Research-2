"""Enclose known compulsory lenses in a common rational polygon at a fixed angle."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import cos,sin,pi
from pathlib import Path
import sys,json
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_compulsory_C3 import EXTRA

def clip(p,n,b):
    h=lambda v:n[0]*v[0]+n[1]*v[1]-b
    res=[]
    for x,y in zip(p,p[1:]+p[:1]):
        a,c=h(x),h(y)
        if a<=0:res.append(x)
        if a*c<0:
            res.append(tuple((a*w-c*v)/(a-c) for v,w in zip(x,y)))
    return res

def run(comp,left,right):
    sys.path.insert(0,str(comp))
    import audit_optimized_eh_chord as g
    from certify_sigma_mandatory_arcs import degree,anchors
    from certify_three_sigma_area_budget import exact,pack
    ctx.prec=192;g.PI=arb.pi();g.SIGMA=degree(left).union(degree(right))
    n,m,*_=g.setup()
    cc=g.line_pair(n[0],m[5]);ll=g.minus(cc,m[5])
    v=g.plus(g.scale(g.HALF,n[1]),g.scale(exact('-.04'),g.rot(n[1])))
    w=g.minus(v,n[1]);j=lambda p:(-p[0],p[1])
    pp=[cc,ll,v,w]+[tuple(exact(c) for c in p) for p in EXTRA]
    defs={k:list(anchors(k)[0].values()) for k in ('TI','GK','QY')}
    defs['C3']=pp+list(map(j,pp))+[(arb(0),exact('.501'))]
    defs['Reuleaux']=[g.line_pair(n[i],n[j]) for i,j in ((0,1),(2,3),(4,5))]
    # Center averaging of the two squared disk inequalities gives an upper
    # support bound for every choice of pair and every rational weight.
    designs={}
    for k,cc0 in defs.items():
        pairs=list(combinations(range(len(cc0)),2))
        lefts=np.array([[float(v) for v in cc0[i]] for i,j in pairs])
        rights=np.array([[float(v) for v in cc0[j]] for i,j in pairs])
        diff=rights-lefts
        d2=np.sum(diff*diff,axis=1)
        designs[k]=(pairs,lefts,diff,d2)
    # The integer directions are a design choice. All actual geometry and all clipping below are exact or directed.
    aa=[360*i/360 for i in range(360)]+[j*60+30+float((left+right)/2) for j in range(6)]
    bb=10**9
    dirs=sorted(set((round(bb*cos(t*pi/180)),round(bb*sin(t*pi/180))) for t in aa))
    polys=[(F(-1),F(-1)),(F(1),F(-1)),(F(1),F(1)),(F(-1),F(1))]
    lines=[]
    for normal in dirs:
        length=arb(normal[0]**2+normal[1]**2).sqrt()
        u=tuple(arb(v)/length for v in normal)
        maxes=[length/2]
        for label,centers in defs.items():
            pairs,lefts,diff,d2=designs[label]
            # Optimize only to choose a rigorous test, never to certify the bound.
            uv=np.array([float(z) for z in u])
            v=diff@uv
            base=lefts@uv
            # Solve the derivative equation for the minimum of
            # dot(u,a+t(b-a)) + sqrt(1-t(1-t)|a-b|^2).
            centered=-v*np.sqrt((1-d2/4)/(d2*np.maximum(d2-v*v,1e-30)))
            weights=np.clip(.5+centered,0,1)
            values=base+v*weights+np.sqrt(1-d2*weights*(1-weights))
            ix=int(np.argmin(values));i,j=pairs[ix]
            t=exact(F(round(float(weights[ix])*10**7),10**7))
            a,b=centers[i],centers[j]
            d=g.minus(b,a);ctr=g.plus(a,g.scale(t,d))
            radius=(1-t*(1-t)*g.dot(d,d)).sqrt()
            maxes.append((g.dot(normal,ctr)+length*radius).upper())
        off=int(max(x.upper() for x in maxes).ceil().unique_fmpz())
        polys=clip(polys,normal,off)
        lines.append([*normal,off])
    area=sum(p[0]*q[1]-p[1]*q[0] for p,q in zip(polys,polys[1:]+polys[:1]))/2
    lim=F('.844082046521218581802948589770767704354')
    return dict(interval=[str(left),str(right)],precision=ctx.prec,outer_polygon_area=pack((area,area),44),
      rigid_disk_intersections=list(defs),num_vertices=len(polys),
      below_B=area<lim,difference_to_B=pack((lim-area,lim-area),36),
      verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('comp',type=Path);p.add_argument('left',type=F);p.add_argument('right',type=F);p.add_argument('slabs',type=int);p.add_argument('out',type=Path)
    a=p.parse_args();rows=[]
    for i in range(a.slabs):
        ends=[a.left+(a.right-a.left)*j/a.slabs for j in (i,i+1)]
        r=run(a.comp.resolve(),*ends)
        rows.append(dict(interval=r['interval'],outer_area=r['outer_polygon_area']['decimal'],below_B=r['below_B'],vertices=r['num_vertices']))
    maximum=max(F(q['outer_area'][1]) for q in rows)
    results=dict(interval=[str(a.left),str(a.right)],all_below_B=all(q['below_B'] for q in rows),
      maximum_upper=rows[max(range(len(rows)),key=lambda i:F(rows[i]['outer_area'][1]))]['outer_area'][1],slabs=rows,
      verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    a.out.write_text(json.dumps(results,indent=1)+'\n');print(results['all_below_B'],results['maximum_upper'])
