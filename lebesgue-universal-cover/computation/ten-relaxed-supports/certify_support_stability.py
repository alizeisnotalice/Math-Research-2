"""Quantify the simultaneous relaxation of six hexagon and four Gibbs supports."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys
from flint import arb,arb_series,ctx,__version__ as flint_version

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def run(comp,prior,tilt):
    sys.path[:0]=[str(comp.resolve()),str(prior.resolve()),str(tilt.resolve())]
    from certify_low_sigma_ah_budget import model,exact,b,g,add,sub,mul,dot,rot,OTHER,pull
    from certify_sigma_mandatory_arcs import anchors,degree
    from certify_sigma_paired_rigidity import rig
    from certify_three_sigma_area_budget import pack
    ctx.prec=192;g.PI=arb.pi()
    W=exact('.00000965');T=exact('.000010')
    eps=exact('.000000011')
    slant=exact('.00000027');line=eps;tangent=eps
    weight0=exact('.000045');tau0=exact('.00004')
    priorwall=tilt/'sharper-wall.json';J=json.loads(priorwall.read_text())
    receipt=comp/'two-moving-boundaries/continuous-joint.json'
    area_receipt=prior/'H-corner-gain.json'
    assert J['sources'][receipt.name]==digest(receipt) and J['sources'][area_receipt.name]==digest(area_receipt)
    source=inspect.getsource(rig)
    marker='    assert min(der)>0'
    assert source.count(marker)==1
    ns=rig.__globals__.copy()
    exec(compile(source[:source.index(marker)]+'    return min(der)\n','endpoint_derivative_projection', 'exec'),ns)
    end=ns['rig']
    lo=F('1.50');hi=F('1.54')
    ws=[F('.435')]+[F(i,20) for i in range(9,21)]
    d=F('.00021');d0=F('.30')
    half=arb(1)/2; rt=arb(3).sqrt(); R=1/rt
    A=exact(F(json.loads(area_receipt.read_text())['full_contour_area_minimum']['exact'][0]))
    floor=lambda key:exact(F(J['lower_bounds'][key]))
    gamma=exact('.0222')
    alpha=3*eps/gamma
    tmax=(2/rt)*(eps+(1-alpha.cos())/2+alpha.sin()/(2*rt))
    drift=(2*eps+alpha/2)/rt
    assert (1-W)/(1-weight0)*(3*weight0/16)-exact('.24')*W>3*eps
    assert (1-T)/(1-tau0)*exact('.24')*tau0-exact('.24')*T>3*eps
    assert (1-T)*exact(F(J['prior_tangent_and_pose_floor']))-3*T > slant+drift
    S={}
    def rec(key,val):
        S[key]=min(S.get(key,val.lower()),val.lower())
    tags=('TI','GK','QY','C3')
    tiny=degree(d);one=degree(F(1))
    for label in tags:
        old=floor(label+'_end')
        angular=(old+9*weight0/16)/(1-weight0)
        B=(angular+3*half*one.sin())/one.cos()
        for angle,needed in ((one,3*eps/tiny),(tiny,gamma)):
            slope=(1-W)*(B*angle.cos()-3*half*angle.sin())-rt*W
            assert slope>needed
            rec('four_tiny' if angle is tiny else 'four_large',slope)
        oldcap=floor(label+'_strong_slant')
        assert (1-W)*(oldcap+weight0/16)/(1-weight0)-W*(R-half)>slant+drift
    # The centered down-pointing Reuleaux triangle has vertices 30,150,270 degrees.
    tri=[(half,1/(2*rt)),(-half,1/(2*rt)),(arb(0),-R)]
    def tri_height(u):
        # It suffices to use a vertex, even where the curved side is active.
        return max(dot(u,p).lower() for p in tri)
    cap_room=W*(degree(hi).cos()/rt-half)
    assert cap_room>slant+drift
    rec('weak_slant',cap_room-drift)
    # Every integer-sector orientation of the prime triangular carrier is pinned
    # to alpha: on [0,pi/6], |F|=sqrt(3)cos(pi/6-t)-3/2 >= .4 t.
    assert (rt-3*half)/(g.PI/6)>exact('.4') and exact('.4')>gamma
    # The AH positive compatibility sum of the triangle on [.1 deg,59.9 deg].
    ah_inside=W*(rt*(g.PI/6-degree(d0)).cos()-3*half)
    assert ah_inside>3*eps
    raw_pivot=floor('paired_tangent_per_triangle_weight')
    assert T*raw_pivot>0
    rows=[]
    for j in range(200):
        l=lo+(hi-lo)*j/200;r=lo+(hi-lo)*(j+1)/200
        g.SIGMA=degree(l).union(degree(r));n,m,*_=g.setup()
        de=model.D_of(g.SIGMA);ta=add(b.z,de)
        E3=g.line_pair(n[4],m[3]);nu=sub(ta,(-E3[0],E3[1]))
        L=dot(nu,ta)
        # Under the QA reflection of the down triangle, (0,1/sqrt(3)) is a point.
        v=(arb(0),R)
        release=dot(nu,v)-L
        # A support point at the canonical hex normal n_j gives a two-sided
        # pin for n_j.t to accuracy eps+(1-cos alpha)/2.
        def projected_coefficient(u,p,z1,z2):
            aa=g.cross(u,n[2])/g.cross(n[1],n[2]);bb=g.cross(n[1],u)/g.cross(n[1],n[2])
            assert min(aa,bb)>0
            return dot(u,rot(p))-aa*dot(n[1],rot(z1))-bb*dot(n[2],rot(z2))
        def pivot_coefficient(a,z):
            s=degree(a).union(degree(z));p=b.anchors(model.D_of(s),exact(ws[0]))
            qj=b.anchors(model.D_of(s),arb_series([exact(ws[0]),1],2))[6]
            vel=tuple(k[1] for k in qj)
            u=mul(-1/g.norm(vel),rot(vel))
            try:
                v0=abs(projected_coefficient(u,p[6],p[5],p[1])).upper()
                if v0<exact('.042'):return v0
            except AssertionError:
                pass
            assert z-a>F(1,10**8),(a,z)
            mid=(a+z)/2;return max(pivot_coefficient(a,mid),pivot_coefficient(mid,z))
        # The carrier's own three terms have absolute sum <= R+1/3 < 1.
        # All six triangular contacts have tangent coordinates <=1/(2sqrt3).
        cpp=(1-T)*pivot_coefficient(l,r)+T
        def projected_loss(C,angle):
            return C*angle.sin()+(2/rt)*(eps+(1-angle.cos())/2)+R*(1-angle.cos())
        pivot_room=T*raw_pivot-projected_loss(cpp,alpha)
        assert pivot_room>tangent,(l,pivot_room)
        rec('pivot_room',pivot_room)
        for w in ws:
            first=end(l,r,w,w,endcap=d0)
            narrow=end(l,r,w,w,endcap=d)
            assert (1-W)*first>3*eps/tiny and (1-W)*narrow>gamma
            rec('AH_narrow',(1-W)*narrow)
            if w==ws[0]:
                assert (1-T)*first-rt*T>3*eps/tiny
                assert (1-T)*narrow-rt*T>gamma
                rec('pivot_narrow',(1-T)*narrow-rt*T)
            g.SIGMA=degree(l).union(degree(r))
            pts=b.anchors(de,exact(w))
            xx=add(b.z,mul(exact(w),de))
            # QA pulls n1 and n2 back to n1 and n0 of the eight-anchor body.
            z1=pull(pts[5],2,-1);z2=pull(pts[3],2,-1)
            coeff=(1-W)*abs(projected_coefficient(nu,xx,z1,z2)).upper()+W
            angle_w=3*eps/((1-W)*narrow)
            line_room=W*release-projected_loss(coeff,angle_w)
            rec('AH_line_coefficient_negative',-coeff)
            assert line_room>line,(l,r,w,line_room)
            rec('AH_line',line_room)
            # Where vertex support is insufficient, central inradius (1-1/sqrt3)
            # is a second uniform triangular lower bound.
            for turn,hand in OTHER:
                rooms=[]
                for jj in (3,5):
                    u=pull(m[jj],turn,hand)
                    h=max((dot(u,z)-half).lower() for z in pts)
                    if hand==1 and (jj-turn)%6==2:h=max(arb(0),h)
                    hu=max(tri_height(u),1-R)
                    rooms.append(((1-W)*h+W*(hu-half)).lower())
                score=max(rooms)-drift
                assert score>slant,(l,w,turn,hand,score)
                rec('AH_ten',score)
        rows.append([str(l),str(r)])
    # Direct mixed-area integration with the actual triangular curvature, 0 or 1.
    # Three active 60-degree arcs; no atoms since the support is C^1.
    mixed=[]
    for j in range(20):
        l=lo+(hi-lo)*j/20;r=lo+(hi-lo)*(j+1)/20
        g.SIGMA=degree(l).union(degree(r));n,m,*_=g.setup()
        ti=anchors('TI')[0];gk=anchors('GK')[0];qy=anchors('QY')[0]
        c=g.line_pair(n[0],m[5])
        pp=tri+[c,(-c[0],c[1]),qy['C'],gk['G'],ti['T'],ti['E'],qy['Q'],
            sub(gk['G'],n[5]),gk['K'],qy['Y'],add(b.z,model.D_of(g.SIGMA))]
        z=ti['B'];pp.append(sub(z,mul(dot(z,m[3])-half,m[3])))
        num=1024;total=arb(0)
        for sector in (1,3,5):
            for k in range(num):
                a=exact(sector+k/F(num))*g.PI/3
                z=exact(sector+(k+1)/F(num))*g.PI/3
                u=g.unit(a.union(z))
                total+=max(half,max(dot(u,p).lower() for p in pp))
        mixed.append((g.PI*total/(6*num)).lower())
    V=min(mixed)
    area0=(1-W)**2*(A+T*floor('straight_point_area_loss_per_triangle_weight_negative'))+2*W*(1-W)*V+W**2*(g.PI-rt)/2
    core_drift=drift
    # Steiner formula for a containing closed convex C dilated by core_drift.
    perimeter=2*rt*(1+2*eps)
    area=area0-perimeter*core_drift-g.PI*core_drift**2
    B=exact('.844082046521218581802948589770767704354')
    assert area>B,(area,area-B)
    for k,z in [('triangle_mixed_area',V),('area_floor',area),('gap_B',area-B),('unmoved_area',area0),
                ('negative_translation_bound',-tmax),('negative_rotation_bound',-alpha),('negative_core_error',-core_drift)]:
        rec(k,z)
    subeps=exact('.000000010');padding=exact('.000000001')
    grid=exact('.00000000001')
    assert F('.000000010')+F('.000000001')==F('.000000011')
    inflate_cost=2*rt*(1+2*subeps)*padding+g.PI*padding**2
    assert area-B>inflate_cost
    rec('finite_net_area_excess_after_inflation',area-B-inflate_cost)
    rec('finite_net_clearance',padding-5*grid)
    return dict(complete=True,interval=['1.50','1.54'],relaxed_distances=['1.1e-8']*6+['2.7e-7']*2+['1.1e-8','1.1e-8'],
                triangle_weight=str(F('.00000965')),pivot_weight=str(F('.000010')),pose_slope='.0222',angular_slices=rows,
                endpoint_bounds=2*len(rows)*len(ws),mixed_direction_boxes=20*3*1024,lower_bounds={k:pack(x)['exact'][0] for k,x in S.items()},
                sources={p.name:digest(p) for p in (priorwall,receipt,area_receipt,comp/'certify_sigma_paired_rigidity.py',Path(__file__))},
                arb_bits=ctx.prec,python=sys.version,python_flint=flint_version)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for key in ('comp','prior','tilt','output'):ap.add_argument(key,type=Path)
    a=ap.parse_args();data=run(a.comp,a.prior,a.tilt);a.output.write_text(json.dumps(data,indent=1)+'\n')
    print({k:float(F(v)) for k,v in data['lower_bounds'].items()})
