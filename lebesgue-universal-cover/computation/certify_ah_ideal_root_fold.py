"""Locate the nontrivial ideal-hexagon compatibility fold and certify the number of interior roots on both sides."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import ctx
from certify_ah_two_hit_budgets import geometry


def run(archive,out):
    ns,_,_=geometry(archive,[(F('.427'),F('.0048'))],far=F('.183'),negative=F('.210'),eps=F('.0051'))
    A=ns['A'];dot=ns['dot'];perp=ns['perp'];cross=ns['cross'];trig=ns['trig'];rot=ns['rotate'];box=ns['box'];normals=list(map(ns['av'],ns['d']['n']));n0=normals[0]
    lo,hi=F('.428'),F('.430');p,q=F('5.5'),F(7)
    selected=[(2,0),(3,1),(0,1)]
    def pts(a,b):
        v=ns['centered_points'](a,b);eds=[ns['sub'](v[k],v[l]) for k,l in ns['orders']]
        return v,eds
    def coeff(v):
        aa=sum((-1)**k*dot(v[i],normals[k]) for k,(i,c) in enumerate(selected))
        bb=sum((-1)**k*dot(v[i],perp(normals[k])) for k,(i,c) in enumerate(selected))
        return aa,bb
    Qmin=A(5);Qmax=-Qmin;unif=A(1);outside=0;peak_bounds=[]
    for wi in range(50):
        a=lo+(hi-lo)*F(wi,50);b=a+(hi-lo)/50;v,e=pts(a,b);w=box(a,b)
        assert all(cross(e[k],e[(k+1)%10])>0 and dot(e[k],e[(k+1)%10])>A('.75') for k in range(10))
        for k,item in enumerate(selected):
            u=rot(n0,trig(box(p+60*k,q+60*k)));cr=[cross(x,u) for x in e];possible=[z for m,z in enumerate(ns['sector']) if not cr[m]<0 and not cr[(m+1)%10]>0]
            assert possible==[item],(a,k,possible)
        aa,bb=coeff(v);da,db=coeff(ns['point_jets'](w)[1]);dq=2*(aa*da+bb*db)
        assert aa>0 and bb>0 and dq<0
        Qmin=min(Qmin,dq.lower());Qmax=max(Qmax,dq.upper());peak=(bb/aa).atan()/ns['angle_unit']
        assert peak>A(p) and peak<A(q);peak_bounds.append(peak)
        slope=sum((-1)**k*dot(v[i],perp(normals[k])) for k,i in enumerate((1,3,0)))
        # From 0 to .5 degrees and symmetrically from 59.5 to 60, F<0.
        assert slope+A('1.5')*(A('.25')*ns['angle_unit']).tan()<0
        def outside_cell(l,r,lev):
            nonlocal outside,unif
            m=(l+r)/2;ss=trig(box(l,r));fm=ns['F_at'](v,e,m)
            dw=ns['angle_unit']*box(l-m,r-m)
            derivative=sum((-1)**k*ns['supports'](v,e,trig(box(l+60*k,r+60*k)),n0,True) for k in range(3))
            gap=(-fm-dw*derivative).lower()
            if gap<=0:
                assert lev<10,(a,b,l,r,gap);outside_cell(l,m,lev+1);outside_cell(m,r,lev+1)
            else:
                outside+=1;unif=min(unif,gap)
        for l,r in ((F('.5'),p),(q,F('59.5'))):
            size=int((r-l)*4)
            for t in range(size):outside_cell(l+(r-l)*F(t,size),l+(r-l)*F(t+1,size),0)
    def discr(w):
        a,b=coeff(pts(w,w)[0]);return a*a+b*b-A('.25')
    assert discr(lo)>0 and discr(hi)<0
    l,r=lo,hi
    for _ in range(88):
        m=(l+r)/2;val=discr(m)
        if val>0:l=m
        elif val<0:r=m
        else:raise ArithmeticError((m,val))
    a,b=coeff(pts(l,r)[0]);fold_angle=(b/a).atan()/ns['angle_unit']
    doc={'w_window':[str(lo),str(hi)],'unique_fold_exact_enclosure':[str(l),str(r)],'fold_w_interval':str(box(l,r)),'fold_angle_degrees':str(fold_angle),'angular_chamber_degrees':[str(p),str(q)],'sector_vertices_and_arc_flags':selected,'discriminant_derivative':[str(Qmin),str(Qmax)],'discriminant_endpoints':[str(discr(lo)),str(discr(hi))],'cell_count_excluding_endpoint_neighborhoods_and_fold_chamber':outside,'strict_negativity_margin_in_cells':str(unif),'certified_angular_end_neighborhoods_degrees':['(0,.5]','[59.5,60)'],'bits':ctx.prec,'python':sys.version,'flint':ns['flint_version'],'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(doc,indent=2)+'\n');print(doc['fold_w_interval'],doc['fold_angle_degrees'],outside)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.output)
