"""Directed boxes for two movable Gibbs slants and two movable AH supports."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
from flint import arb, arb_series, ctx, __version__ as flint_version

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def run(comp, prior):
    sys.path.insert(0,str(Path(__file__).parent))
    from wall_refinement import strengthened_root, fourth, mixed_area
    refinement=Path(__file__).parent/'wall_refinement.py'
    sys.path[:0]=[str(comp.resolve()),str(prior.resolve())]
    from certify_low_sigma_ah_budget import model,exact,b,g,add,sub,mul,dot,rot,OTHER,pull
    from certify_sigma_mandatory_arcs import degree,anchors
    from certify_compulsory_C3 import EXTRA
    from certify_three_sigma_area_budget import pack
    ctx.prec=192
    g.PI=arb.pi()
    half=arb(1)/2
    weight=exact('.000045')
    epsilon=exact('.0000028')
    eta=exact('.000000123')
    tau=exact(F(4,100000))
    receipt=comp/'two-moving-boundaries/continuous-joint.json'
    area_receipt=prior/'H-corner-gain.json'
    joint=json.loads(receipt.read_text())
    left=F('1.50');right=F('1.54')
    assert joint['complete'] and joint['interval']==[str(left),str(right)]
    root_min=min(F(r['root_margin_lower']['exact'][0]) for r in joint['panels_data'])
    pose_min=min(F(r['discrete_pose_lower']['exact'][0]) for r in joint['panels_data'])
    # The stronger pivot-only root check in fourth() carries the triangle drift.
    assert (1-tau)*exact(pose_min)-3*tau>epsilon
    rerig=strengthened_root()
    mixed_values=[mixed_area(left+(right-left)*j/20,left+(right-left)*(j+1)/20) for j in range(20)]
    V=min(mixed_values)
    ws=[F('.435')]+[F(i,20) for i in range(9,21)]
    panes=[];stats={}
    def record(k,v):stats[k]=min(stats.get(k,v.lower()),v.lower())
    def harmonic(u):return (3*u[1]-4*u[1]**3)/16
    for i,row in enumerate(joint['panels_data']):
        l,r=map(F,row['degrees'])
        assert (l,r)==(left+(right-left)*i/200,left+(right-left)*(i+1)/200)
        g.SIGMA=degree(l).union(degree(r));n,m,*_=g.setup()
        root_count,cost,tangent=fourth(l,r,g.SIGMA,rerig)
        record('straight_point_area_loss_per_triangle_weight_negative',-cost)
        record('paired_tangent_per_triangle_weight',tangent)
        weak=weight*(3*g.SIGMA).cos()/16
        assert weak>epsilon;record('weak_slant',weak)
        D=model.D_of(g.SIGMA)
        def reflected_gap(s):
            ta=add(b.z,model.D_of(s))
            e3=g.line_pair(n[4],g.unit(7*g.PI/6+s));out=sub(ta,(-e3[0],e3[1]))
            return half+harmonic(pull(out,2,-1))-dot(out,ta)
        middle=degree((l+r)/2)
        gap=reflected_gap(middle)+reflected_gap(arb_series([g.SIGMA,1],2))[1]*(g.SIGMA-middle)
        assert weight*gap>eta,(l,r,gap);record('AH_line_excess',weight*gap)
        for w in ws:
            points=b.anchors(D,exact(w))
            for k,hand in OTHER:
                excess=[]
                for j in (3,5):
                    u=pull(m[j],k,hand)
                    floor=max((dot(u,v)-half).lower() for v in points)
                    if hand==1 and (j-k)%6==2:floor=max(floor,arb(0))
                    if floor>=0:
                        excess.append(((1-weight)*floor+weight*harmonic(u)).lower())
                best=max(excess)
                assert best>epsilon, (l,w,k,hand,best)
                record('AH_ten_slants',best)
        c=g.line_pair(n[0],m[5])
        v=add(mul(half,n[1]),mul(exact('-.04'),rot(n[1])))
        cp={'C':c,'L':sub(c,m[5]),'V':v,'W':sub(v,n[1]),'A':(arb(0),arb('.501'))}
        cp.update({f'k{j}':tuple(exact(z) for z in p) for j,p in enumerate(EXTRA)})
        cp.update({'J'+k:(-p[0],p[1]) for k,p in list(cp.items()) if k!='A'})
        counts={'pivot_triangle':root_count}
        for label in ('TI','GK','QY','C3'):
            if label=='C3':
                points,phase,trip,guard=cp,1,('JC','JW','V'),(0,1)
            else:
                points,_,phase,trip,raw=anchors(label)
                guard=tuple(j for j,_ in raw)
            # At the even remaining weak normals a unit diameter contact
            # gives support >=1/2. All other positive normals are guarded below.
            assert all(j%2==0 for j in set((0,1,2)).difference(guard))
            cap=min(((1-weight)*max((dot(m[j],p)-half).lower() for p in points.values())-weight/16).lower() for j in guard)
            assert cap>epsilon
            record(label+'_strong_slant',cap)
            def dirs(l,r):
                arg=degree(l).union(degree(r))
                return [g.unit(phase*g.PI+arg+2*j*g.PI/3) for j in range(3)]
            slope0=sum(dot(rot(u),points[k]) for u,k in zip(dirs(F(0),F(0)),trip))
            deriv=(1-weight)*(slope0*degree(F(1)).cos()-3*half*degree(F(1)).sin())-9*weight/16
            assert deriv>0
            record(label+'_end',deriv)
            todo=[(F(1),F(30))];num=0
            while todo:
                a,z=todo.pop()
                floor=sum(max(dot(u,p).lower() for p in points.values()) for u in dirs(a,z))-3*half
                room=(1-weight)*floor-3*weight/16
                if room>0:
                    record(label+'_root',room);num+=1;continue
                assert z-a>F(1,10**7),(label,l,r,a,z,room)
                mm=(a+z)/2;todo.extend(((a,mm),(mm,z)))
            counts[label]=num
        panes.append(dict(sigma=row['degrees'],root_cells=counts))
    assert len(panes)==200
    ar=json.loads(area_receipt.read_text());assert ar['complete'] and ar['interval']==joint['interval']
    assert ar['dependencies']['continuous-joint.json']==digest(receipt)
    A=exact(F(ar['full_contour_area_minimum']['exact'][0]))
    A+=tau*stats['straight_point_area_loss_per_triangle_weight_negative']
    area=(1-weight)**2*A+2*weight*(1-weight)*V+weight**2*(15*g.PI/64)
    B=exact('.844082046521218581802948589770767704354')
    assert area-B>exact('.0000000114')
    return dict(complete=True,interval=joint['interval'],positive_caps=['.0000028','.0000028'],
                straight_AH='.000000123',paired_AH_tangent='.00000036',weight=str(F('.000045')),
                triangle_weight=str(F(4,100000)),prior_root_floor=str(root_min),prior_tangent_and_pose_floor=str(pose_min),
                lower_bounds={k:pack(v)['exact'][0] for k,v in stats.items()},
                mixed_area_lower=pack(V)['exact'][0],area_lower=pack(area)['exact'][0],excess_over_B=pack(area-B)['exact'][0],cells=panes,
                sources={p.name:digest(p) for p in (receipt,area_receipt,prior/'certify_compulsory_C3.py',
                comp/'certify_sigma_mandatory_arcs.py',comp/'certify_low_sigma_ah_budget.py',comp/'certify_sigma_paired_rigidity.py',
                comp/'two-moving-boundaries/certify_continuous_joint.py',Path(__file__),refinement)},
                python=sys.version,flint=flint_version,precision=ctx.prec)

if __name__=='__main__':
    ap=ArgumentParser()
    for k in ('comp','prior','output'):ap.add_argument(k,type=Path)
    a=ap.parse_args();res=run(a.comp,a.prior)
    a.output.write_text(json.dumps(res,indent=1)+'\n')
    print({k:float(F(v)) for k,v in res['lower_bounds'].items()},float(F(res['area_lower'])),float(F(res['excess_over_B'])))
