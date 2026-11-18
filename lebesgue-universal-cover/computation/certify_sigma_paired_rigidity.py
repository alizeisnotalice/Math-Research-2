"""Rigidity of a continuum of paired AH points with only the old Ta-Z half-plane retained."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

from certify_sigma_mandatory_arcs import degree,g,rotate_angle


def rig(s0,s1,w0,w1,*,endcap=F(1),min_cell=F(1,1024),probe_angles=None):
    h=g.HALF
    g.SIGMA=degree(s0).union(degree(s1))
    n,m,*_=g.setup()
    rt=g.arb(3).sqrt()
    from flint import arb_series
    mid=(-h/2,-rt/4)
    z=g.circle_pair((rt/2-h,(rt/2-h-1)/rt),mid,-1)
    qa=lambda v:rotate_angle((v[0],-v[1]),2*g.PI/3)
    el=(rt/2,h)
    en=(-rt/2,h)
    def raw(sigma,w):
        e3=g.line_pair(n[4],g.unit(7*g.PI/6+sigma))
        c=(-e3[0],e3[1])
        v=g.minus(z,c)
        vv=g.dot(v,v)
        ta=g.plus(c,g.scale(1/vv,g.minus(v,g.scale((vv-1).sqrt(),g.rot(v)))))
        D=g.minus(ta,z)
        out=g.minus(ta,c)
        dy=qa(out)[1]
        removal=h+(3*dy-4*dy**3)/16-g.dot(out,ta)
        x=g.plus(z,g.scale(w,D))
        def inverse(p,normal,direction):
            off=g.dot(p,normal)-h
            return g.minus(g.minus(p,g.scale(off,normal)),g.scale((1-off**2).sqrt(),direction))
        ll=inverse(x,n[5],el)
        nt=inverse(x,n[4],en)
        u=ll[0]+h
        xx=(u+rt*(1-u*u).sqrt())/2-h
        b=(xx,(xx-1)/rt)
        f,N=qa(ll),qa(nt)
        yp=g.circle_pair(b,N,-1)
        ps=[b,g.minus(b,n[5]),f,g.plus(f,n[0]),N,g.minus(N,n[4]),yp,qa(x)]
        return [v for p in ps for v in p]+[removal,g.dot(g.minus(x,c),el)]
    sm=degree((s0+s1)/2)
    wm=g.exact((w0+w1)/2)
    sb=degree(s0).union(degree(s1))
    wb=g.exact(w0).union(g.exact(w1))
    sr=degree((s1-s0)/2)
    wr=g.exact((w1-w0)/2)
    def jet(v,d):
        return arb_series([v,d],2)
    central=raw(sm,wm)
    sder=raw(jet(sb,sr),jet(wb,g.arb(0)))
    wder=raw(jet(sb,g.arb(0)),jet(wb,wr))
    combined=[]
    for center,ds,dw in zip(central,sder,wder,strict=True):
        assert ds[1].is_finite() and dw[1].is_finite()
        combined.append(center+ds[1].union(-ds[1])+dw[1].union(-dw[1]))
    ps=list(zip(combined[0:16:2],combined[1:16:2],strict=True))
    removal,branch=combined[16:]
    assert min(removal,branch)>0
    b,yp=ps[0],ps[6]
    if w1==1:
        vr=g.minus(yp,b)
        assert g.cross(vr,n[2])>0 and g.dot(vr,n[2])>0
    units={(0,1),(2,3),(4,5),(0,2),(0,6),(2,7),(4,6),(4,7)}
    slack=min((1-g.dot(g.minus(ps[i],ps[j]),g.minus(ps[i],ps[j]))).lower()
              for i in range(8) for j in range(i+1,8) if (i,j) not in units)
    assert slack>0
    mp=[g.unit(g.PI/6-g.SIGMA+i*g.PI/3) for i in range(6)]
    caps=[g.dot(m[0],ps[7])-h,g.dot(m[1],yp)-h,g.dot(mp[0],ps[7])-h,
          g.dot(mp[1],yp)-h,g.dot(mp[5],b)-h]
    assert min(caps)>0 and g.dot(mp[5],el)>0 and (3*g.SIGMA).cos()>0
    verts=[ps[i] for i in (0,2,7,4,6)]
    wedges=[]
    for i,p in enumerate(verts):
        a=g.minus(p,verts[(i-1)%5])
        b0=g.minus(p,verts[(i+1)%5])
        if g.cross(a,b0)<0:a,b0=b0,a
        assert g.cross(a,b0)>0
        wedges.append((a,b0,p))
    @lru_cache(maxsize=None)
    def directions(a,b):
        angle=degree(a).union(degree(b))
        return [g.unit(angle+i*2*g.PI/3) for i in range(3)]
    def support(u):
        ans=max(g.dot(u,p).lower() for p in ps)
        d=g.neg(u)
        for a,b,p in wedges:
            if g.cross(a,d)>=0 and g.cross(d,b)>=0:
                ans=max(ans,(1+g.dot(u,p)).lower())
        return ans
    if probe_angles is not None:
        return [(sum(support(v) for v in directions(a,a))-3*h).lower() for a in probe_angles]
    der=[]
    for theta,sign,ix in [(0,1,[3,1,4]),(60,-1,[5,2,0])]:
        B=sum(sign*g.dot(g.rot(v),ps[j]) for j,v in zip(ix,directions(F(theta),F(theta)),strict=True))
        der.append((B*degree(endcap).cos()-3*h*degree(endcap).sin()).lower())
    assert min(der)>0
    vals=[]
    def scan(a,b):
        gap=sum(support(v) for v in directions(a,b))-3*h
        if gap>0:
            vals.append(gap.lower())
            return
        assert b-a>min_cell
        k=(a+b)/2
        scan(a,k)
        scan(k,b)
    scan(endcap,F(60)-endcap)
    return len(vals),[t.lower() for t in (slack,min(caps),removal,min(der),min(vals))]


def verify(low,high,wmin,wmax,ns,nw):
    assert F(1)<=low<high<=F(3) and 0<=wmin<=wmax<=1 and min(ns,nw)>0
    records=[]
    for j in range(ns):
        a,b=[low+(high-low)*i/ns for i in (j,j+1)]
        for k in range(nw):
            c,d=[wmin+(wmax-wmin)*i/nw for i in (k,k+1)]
            n,res=rig(a,b,c,d)
            records.append((n,res))
    return {'sigma_degrees':[str(low),str(high)],'paired_parameter':[str(wmin),str(wmax)],'sigma_cells':ns,
            'paired_cells':nw,'joint_root_cells':sum(q for q,_ in records),
            'lower_bounds':{k:str(min(row[i] for _,row in records)) for i,k in enumerate(
                ['nonunit_diameter_slack','five_strict_caps','reflected_harmonic_excess','endpoint_root_derivative','interior_root_excess'])},
            'precision_bits':g.ctx.prec}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('low','high','wmin','wmax'):p.add_argument(key,type=F)
    for key in ('ns','nw'):p.add_argument(key,type=int)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(verify(args.low,args.high,args.wmin,args.wmax,args.ns,args.nw),indent=2)+'\n')
