"""Check whether the entire 21-anchor lens adds anything to the H contour."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sys,json
from flint import arb,arb_series,ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_compulsory_C3 import EXTRA

def run(comp):
    sys.path.insert(0,str(comp))
    from certify_low_sigma_ah_budget import g,b,model,add,sub,mul,cross,dot,rot
    from certify_three_sigma_area_budget import exact,pack
    ctx.prec=192
    mins=[];cuts=[];slabs=200
    def visit(left,right):
        s=b.degree(left).union(b.degree(right))
        n=b.directions;m=[g.unit(g.PI/6+s+i*g.PI/3) for i in range(6)]
        c=g.line_pair(n[0],m[5]);ell=sub(c,m[5])
        v=add(mul(g.HALF,n[1]),mul(exact('-.04'),rot(n[1])))
        w=sub(v,n[1]);j=lambda p:(-p[0],p[1])
        x=[c,ell,v,w]+[tuple(exact(v0) for v0 in p) for p in EXTRA]
        x+=list(map(j,x))+[(arb(0),exact('.501'))]
        h=j(c);bb=g.line_pair(n[1],m[0]);r=g.line_pair(m[2],n[1],g.HALF,0)
        ss=g.circle_line(r,n[0],-1);t=g.circle_pair(ss,bb,1)
        e=sub(ss,n[0]);a=sub(h,c);bb=sub(h,j(ell));normal=mul(-1,rot(sub(t,h)))
        signs=[cross(a,normal),cross(normal,bb),cross(a,bb),cross(e,h),cross(h,t)]
        assert min(signs)>0
        mid=b.degree((left+right)/2);width=b.degree((right-left)/2)
        D=model.D_of(mid);Dj=model.D_of(arb_series([s,1],2));ws=[F('.435')]+[F(i,20) for i in range(9,21)]
        cp0=add(b.z,mul(exact(ws[0]),D));cpj=add(b.z,mul(exact(ws[0]),Dj))
        def bound(p,pj):
            return tuple(v+(-width).union(width)*germ[1] for v,germ in zip(p,pj,strict=True))
        cp=bound(cp0,cpj)
        qs=[bound(b.anchors(D,exact(w))[6],b.anchors(Dj,exact(w))[6]) for w in ws]
        ds=[]
        for p,q in zip([cp]+qs,qs):
            edge=sub(q,p);u=mul(-1/g.norm(edge),rot(edge))
            # Every anchor yields its disk support upper; pick the best valid disk.
            margin=max((dot(u,sub(p,z))-1).lower() for z in x)
            assert margin>0
            ds.append(margin)
        return min(z.lower() for z in signs),min(ds)
    def test(a,z):
        try: sg,room=visit(a,z)
        except AssertionError:
            assert z-a>F(1,10000000)
            mid=(a+z)/2;test(a,mid);test(mid,z)
        else:mins.append(sg);cuts.append(room)
    for i in range(slabs):
        a=F('1.5')+F(i,5000);test(a,a+F(1,5000))
    return dict(interval=['1.50','1.54'],complete=True,closed_slabs=len(mins),AH_new_supports=13,
        min_cross=pack(min(mins),35),min_AH_support_room=pack(min(cuts),35),
        verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('comp',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();r=run(a.comp.resolve());a.out.write_text(json.dumps(r,indent=1)+'\n')
    print(r['closed_slabs'],r['min_AH_support_room']['decimal'],r['min_cross']['decimal'])
