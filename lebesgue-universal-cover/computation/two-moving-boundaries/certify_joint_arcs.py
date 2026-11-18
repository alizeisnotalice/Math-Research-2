"""Lower bound with freely redesigned, simultaneous AH and EH boundaries inside the original Gibbs host."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json, subprocess, sys

from flint import arb, arb_series, ctx


def run(comp):
    sys.path.insert(0,str(comp))
    import flint
    from certify_low_sigma_ah_budget import model, exact, b, g, add, sub, mul, dot, cross, rot, neg, OTHER, pull
    from certify_three_sigma_area_budget import setup, signed_arc, pack
    from certify_sigma_mandatory_arcs import anchors, certify_cell
    from certify_sigma_paired_rigidity import rig

    ctx.prec=192
    ctx.cap=10
    threshold=F('.844082046521218581802948589770767704354')
    pivot=F('.435')
    nodes=[pivot]+[F(i,20) for i in range(9,21)]
    results=[]
    for text in ('1.50','1.52','1.53','1.54'):
        s=F(text)
        V=setup(text,g)
        D=model.D_of(b.degree(s))
        straight=add(b.z,mul(exact(pivot),D))
        # Isolate all possible continuous rotations for thirteen actual anchors.
        roots=[]
        for w in nodes:
            n,stats=rig(s,s,w,w,endcap=F('.1'),min_cell=F(1,1000000))
            assert min(stats)>0
            roots.append(dict(w=str(w),closed_angle_boxes=n,lower_bounds=[pack(z) for z in stats]))
        pp=b.anchors(D,exact(pivot))
        x=pp[6]
        curve=b.anchors(D,arb_series([exact(pivot),1],2))[6]
        velocity=tuple(z[1] for z in curve)
        nu=mul(-1/g.norm(velocity),rot(velocity))
        inward=rot(mul(-1/g.norm(D),D))
        vert=(arb(0),1/arb(3).sqrt())
        margins=[dot(nu,sub(vert,x)),-dot(inward,sub(vert,straight))]
        n,m,*_=g.setup()
        for turn,hand in OTHER:
            outward=[pull(m[i],turn,hand) for i in (3,5)]
            margins.append(max((dot(v,a)-g.HALF).lower() for v in outward for a in pp))
        assert len(margins)==12 and min(margins)>0
        # Signed Green area between old AH and the chain of compulsory points.
        def rel(w): return sub(b.anchors(D,w)[6],straight)
        def density(t):
            xx,yy=rel(arb_series([t,1],8))
            return xx*yy.derivative()-yy*xx.derivative()
        acc=arb(0)
        cells=192
        h=exact(F(1,cells))
        for i in range(cells):
            l=F(i,cells);r=F(i+1,cells)
            coeff=density(exact((l+r)/2));remainder=density(b.ea(l,r))[6]
            acc+=sum(coeff[k]*h**(k+1)/((k+1)*2**k) for k in (0,2,4))+remainder*h**7/(7*2**6)
        poly=sum(cross(rel(exact(a)),rel(exact(z))) for a,z in zip(nodes,nodes[1:]))
        ah=(acc-poly)/2
        assert ah>0

        # The jointly compulsory TI lens extends all the way up to E=S-(1,0).
        ncell,ti=certify_cell('TI',s,s)
        pts=anchors('TI')[0]
        E,S,T,B=(pts[k] for k in ('E','S','T','B'))
        H,U,C3=(V[k] for k in ('H','U','C3'))
        va,vb=sub(E,S),sub(T,S)
        assert cross(va,vb)>0 and dot(va,vb)>0
        angle=arb.atan2(cross(va,vb),dot(va,vb))
        assert 0<angle<exact(F('.04'))
        # An acute unit arc has distance to T at most its angular length.
        c01,s01=g.unit(arb(0).union(angle))
        radial=add(mul(c01,va),mul(s01,rot(va)))
        dB=2*dot(sub(S,B),rot(radial))
        assert dB>0 and ti['lens_distance_slack']>0
        HtoU=sub(H,C3),sub(U,C3)
        assert cross(*HtoU)>0 and dot(*HtoU)>0
        assert E[1]>H[1]>U[1]>T[1]
        # At T the new S-arc meets the old B3-arc; the two tangents turn left.
        assert cross(vb,sub(T,B))>0
        # At E the unit S arc has precisely the downward vertical tangent.
        eh=(cross(E,H)+signed_arc(H,U,C3,g)+cross(U,T)-signed_arc(E,T,S,g))/2
        assert eh>0
        chord=sub(T,E);out=mul(-1/g.norm(chord),rot(chord))
        v_mid=add(va,vb);v_mid=mul(1/g.norm(v_mid),v_mid)
        mid=add(S,v_mid)
        assert cross(va,v_mid)>0 and cross(v_mid,vb)>0
        lam=exact(F(1,10**7))
        excess0=dot(out,sub(mid,E))
        excess=dot(out,sub(mul(1-lam,mid),E))
        assert excess>exact('.00012') and excess0>0
        chord_eh=(cross(E,H)+signed_arc(H,U,C3,g)+cross(U,T)-cross(E,T))/2
        assert chord_eh>eh
        main_area=json.loads(subprocess.check_output([sys.executable,'-B',str(comp/'certify_gibbs_area_refined.py'),text,'16','64']))
        assert main_area['passed'] and F(main_area['sigma_degrees'])==s
        base=tuple(map(F,main_area['area_interval']))
        lower=base[0]-F(pack(ah)['exact'][1])-F(pack(eh)['exact'][1])
        upper=base[1]-F(pack(ah)['exact'][0])-F(pack(eh)['exact'][0])
        assert upper>lower>threshold
        failed_low=base[0]-F(pack(ah)['exact'][1])-F(pack(chord_eh)['exact'][1])
        failed_high=base[1]-F(pack(ah)['exact'][0])-F(pack(chord_eh)['exact'][0])
        assert failed_low<failed_high<threshold
        results.append(dict(sigma=text,rigid_points=roots,positive_pose_excess=[pack(v) for v in margins],
            AH_green_cap_area=pack(ah,32),EH_max_removed=pack(eh,32),gibbs_area=pack(base,32),
            compulsory_enclosed_area=pack((lower,upper),32),gap_above_B=pack((lower-threshold,upper-threshold),32),
            rejected_straight_ET=dict(area=pack((failed_low,failed_high),32),EH_removed=pack(chord_eh,32),
                forced_arc_middle=[pack(v,28) for v in mid],normal=[pack(v,28) for v in out],
                limit_signed_distance=pack(excess0,28),perturbed_point_distance=pack(excess,28),harmonic_weight='1/10000000'),
            TI_interior_angles=ncell,TI_rechecked={k:pack(v) for k,v in ti.items()},
            ET_radians=pack(angle),ET_distance_squared_to_B_derivative_lower=pack(dB.lower())))
    depends=('certify_low_sigma_ah_budget.py','certify_three_sigma_area_budget.py','certify_sigma_paired_rigidity.py',
             'certify_sigma_mandatory_arcs.py','certify_two_sigma_ideal_placement.py','certify_moving_sigma_fold_splice.py',
             'certify_gibbs_area_refined.py','audit_optimized_eh_chord.py')
    return dict(complete=True,target=str(threshold),precision_bits=ctx.prec,flint=flint.__version__,python=sys.version,
                pivot=str(pivot),AH_panels=192,derivative_remainder_order=6,results=results,
                source_sha256={f:sha256((comp/f).read_bytes()).hexdigest() for f in depends},
                verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('computation',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();q=run(a.computation.resolve())
    a.output.write_text(json.dumps(q,indent=2)+'\n')
    for r in q['results']:print(r['sigma'],r['compulsory_enclosed_area']['decimal'],r['gap_above_B']['decimal'])
