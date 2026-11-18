"""Certify exclusion of the Tp-P cut, and necessity of that entire arc when Ta-Z is retained."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_general_relocation import exact_fraction, rotate_angle
from certify_eh_joint_parameters import specialize


def normalize(x):
    return g.scale(1/g.norm(x),x)


def verify():
    specialize(F('1.5497'))
    n,m,*_=g.setup()
    mp=[g.unit(g.PI/6-g.SIGMA+j*g.PI/3) for j in range(6)]
    qa=lambda p:rotate_angle((p[0],-p[1]),2*g.PI/3)
    a=g.arb
    f3=g.line_pair(n[3],m[2])
    e3=g.line_pair(n[4],m[3])
    G=g.circle_line(f3,n[5],1)
    c=(-e3[0],e3[1])
    rt=a(3).sqrt()
    mid=(-g.HALF/2,-rt/4)
    z=g.circle_pair((rt/2-g.HALF,(rt/2-g.HALF-1)/rt),mid,-1)
    v=g.minus(z,c)
    vv=g.dot(v,v)
    ta=g.plus(c,g.scale(1/vv,g.minus(v,g.scale((vv-1).sqrt(),g.rot(v)))))
    segment=g.minus(ta,z)
    outward=normalize(g.rot(segment))
    el=(rt/2,g.HALF)
    en=(-rt/2,g.HALF)
    radial=rotate_angle(n[2],-a('.022'))
    toe=g.plus(G,g.scale(1-a('.000020'),radial))
    q=g.plus(g.plus(z,g.scale(a('.99'),segment)),g.scale(a('.0000002'),outward))
    # c is the vertical mirror of E3, hence qa(c)=F3 exactly by the two defining face lines.
    b=g.minus(G,g.scale(a('.0002'),el))
    f=qa(g.plus(c,g.scale(a('.00008'),el)))
    de=g.plus(mid,g.scale(a('-.096'),en))
    all_points=lambda b,f,de,y,yq,marks:[b,g.minus(b,n[5]),f,g.plus(f,n[0]),de,g.minus(de,n[4]),y,yq,
                                         *(tuple(map(a,row)) for row in marks)]
    strict=all_points(b,f,de,toe,qa(q),[('.10563549','-.46244989'),('-.33973885','.33986545')])
    ideal_de=qa(g.circle_line(ta,n[4],1))
    tp=g.circle_pair(G,ideal_de,-1)
    frozen=all_points(G,f3,ideal_de,tp,qa(ta),[('.10882837','-.45902979'),('-.33956564','.33996545')])
    unit=[(0,1),(2,3),(4,5)]
    extra=[(0,2),(0,6),(2,7),(4,6),(4,7)]
    sec=lambda deg:exact_fraction(deg)*g.PI/180
    ans={}
    for tag,points,ones in [('strict',strict,unit),('limit',frozen,unit+extra)]:
        gaps=[(1-g.dot(g.minus(points[i],points[j]),g.minus(points[i],points[j]))).lower()
              for i,j in combinations(range(len(points)),2) if (i,j) not in ones]
        assert min(gaps)>0
        def dirs(t):
            return [g.unit(t+k*2*g.PI/3) for k in range(3)]
        endpoint=[]
        for angle,side,ids in [(0,1,[3,1,4]),(60,-1,[5,2,0])]:
            # The three designated contacts attain 1/2 on all their directions at this endpoint.
            B=sum(side*g.dot(g.rot(u),points[j]) for u,j in zip(dirs(sec(F(angle))),ids,strict=True))
            bound=B*sec(F(1)).cos()-3*g.HALF*sec(F(1)).sin()
            assert bound>0
            endpoint.append(str(bound.lower()))
        cells=[]
        def check(l,r):
            theta=sec(l).union(sec(r))
            score=sum(max(g.dot(u,v).lower() for v in points) for u in dirs(theta))-3*g.HALF
            if score>0:
                cells.append([str(l),str(r),str(score.lower())])
                return
            assert r-l>F(1,4096)
            h=(l+r)/2
            check(l,h)
            check(h,r)
        check(F(1),F(59))
        guards={}
        for sig,normals,ids in [('+',m,[0,1,2]),('-',mp,[0,1,5])]:
            for i in ids:
                scores=[(g.dot(normals[i],v)-g.HALF).lower() for v in points]
                best=max(scores)
                if tag=='limit' and (sig,i)==('+',2):
                    assert scores.index(best)==2  # defining equality m2.F3=1/2
                else:assert best>0
                guards[f'{sig}{i}']=[scores.index(best),str(best)]
        ans[tag]={'points':[list(map(str,p)) for p in points],'exact_unit_pairs':ones,'other_pairs':len(gaps),
                  'other_squared_distance_slack_lower':str(min(gaps)),'cap_guards':guards,
                  'first_last_degree_derivative_lower':endpoint,'interior_cells':cells,
                  'interior_slack_lower':str(min(g.arb(x[2]).lower() for x in cells))}
    poses=[]
    for k in range(6):
        if not {(3-k)%6,(5-k)%6}&{0,1,2}:poses.append([k,False])
        if not {(k-4)%6,(k-6)%6}&{0,1,5}:poses.append([k,True])
    assert poses==[[0,False],[2,True]]
    new_excess=g.dot(radial,g.minus(toe,G))-a('.0085').cos()
    old_excess=g.dot(outward,g.minus(q,z))
    assert new_excess>0 and old_excess>0
    vr=g.minus(tp,G)
    span=g.arb.atan2(g.cross(vr,n[2]),g.dot(vr,n[2]))
    assert 0<a('.0305')<span<g.PI/2
    return {'sigma_degrees':'1.5497','precision_bits':g.ctx.prec,'cases':ans,'permitted_dihedral_poses':poses,
            'strict_original_new_chord_excess_lower':str(new_excess.lower()),
            'strict_reflected_Ta_Z_excess_lower':str(old_excess.lower()),'whole_Tp_P_sweep_radians':str(span)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(verify(),indent=2)+'\n')
