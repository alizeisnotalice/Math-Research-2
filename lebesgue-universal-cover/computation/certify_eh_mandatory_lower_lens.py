"""A mandatory two-dimensional lower-tip lens for every closed universal subset of P."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_general_relocation import exact_fraction, rotate_angle
from certify_eh_joint_parameters import specialize


def verify():
    specialize(F('1.5497'))
    n,m,_,b,_,_,s,t,_=g.setup()
    q=lambda z:rotate_angle((z[0],-z[1]),5*g.PI/3)
    r=g.line_pair(m[2],n[1],g.HALF,g.arb(0))
    f=g.unit(5*g.PI/6-g.SIGMA)
    a=tuple(map(g.arb,['.125327173657','.467705236099']))
    pts={'S':s,'E':g.minus(s,n[0]),'R':r,'U':g.minus(r,f),'B':b,'Z':g.minus(b,n[1]),'T':t,'a':a}
    for k in ('S','E','U','T','a'):
        pts['Q'+k]=q(pts[k])
    unit_pairs=[('S','E'),('S','R'),('S','T'),('R','U'),('B','Z'),('B','T'),('R','QU'),('R','QS'),('QS','QT'),('QS','QE'),('Z','QT')]
    exact={frozenset(z) for z in unit_pairs}
    gaps=[(1-g.dot(g.minus(v,w),g.minus(v,w))).lower() for (i,v),(j,w) in combinations(pts.items(),2)
          if frozenset((i,j)) not in exact]
    assert min(gaps)>0
    radius=g.arb('.001')
    lens_gap=min((1-g.norm(g.minus(t,z))-radius).lower() for i,z in pts.items() if i not in ('S','B'))
    mirror_gap=1-g.norm(g.minus(t,q(t)))-2*radius
    assert lens_gap>0 and mirror_gap>0

    lam=g.arb('.0000001')  # covers every 0 < mixing coefficient <= lam
    def angles(theta):
        return [theta+g.PI+j*2*g.PI/3 for j in range(3)]
    def witness_bound(theta):
        return sum(max(g.dot(g.unit(p),z).lower() for z in pts.values()) for p in angles(theta))-3*g.HALF
    # At 0, E, QS and B give exactly 1/2 each along n3,n5,n1.
    A=3*g.HALF
    B=sum(g.dot(g.rot(g.unit(p)),pts[k]) for p,k in zip(angles(g.arb(0)),('E','QS','B'),strict=True))
    rad=lambda deg:exact_fraction(deg)*g.PI/180
    angle0=rad(F(1))
    derivative=B*angle0.cos()-A*angle0.sin()
    # Rotated opposite triple has harmonic term -.03 sin(3 theta).
    derivative=(1-lam)*derivative-9*lam/100
    assert derivative>0
    cells=[]
    def partition(start,end):
        theta=rad(start).union(rad(end))
        floor=witness_bound(theta)
        lower=(1-lam)*floor-3*lam/100
        if lower>0:
            cells.append([str(start),str(end),str(lower.lower())])
            return
        assert end-start>F(1,2048),(start,end)
        mid=(start+end)/2
        partition(start,mid)
        partition(mid,end)
    partition(F(1),F(30))

    # Positive cap supports are witnessed symbolically at m0.B=m2.R=1/2;
    # Q supplies the matching prime caps. Middle cap has a strict coordinate margin.
    cap_mid=g.dot(m[1],q(t))-g.HALF
    end_sine=(3*g.SIGMA).cos()/100
    mid_margin=(1-lam)*cap_mid-lam/100
    assert end_sine>0 and mid_margin>0
    # Exhaust D6. Before rotation, even elements use m_(k-shift).
    # Reflected transformations send m_k to m'_(shift-k-1).
    positive=[{0,1,2},{2,3,4}]
    admitted=[]; excluded=[]
    for shift in range(6):
        for refl in (False,True):
            pull=[(shift-k-1 if refl else k-shift)%6 for k in (3,5)]
            hit=sorted(set(pull)&positive[int(refl)])
            if hit:excluded.append({'shift_60deg':shift,'reflection':refl,'cap_preimages':pull,'hit':hit})
            else:admitted.append((shift,refl))
    assert admitted==[(0,False),(5,True)]
    return {'sigma_degrees':'1.5497','precision_bits':g.ctx.prec,'new_rational_anchor':['.125327173657','.467705236099'],
            'anchor_count':len(pts),'exact_unit_pairs':unit_pairs,'other_pair_count':len(gaps),
            'other_squared_distance_slack_lower':str(min(gaps)),'forced_lens_radius':str(radius),
            'other_anchor_margin_after_radius_lower':str(lens_gap),'mirror_lens_diameter_slack_lower':str(mirror_gap.lower()),
            'harmonic_amplitude':'1/100','upper_mixing_weight':'1/10000000','first_degree_derivative_lower':str(derivative.lower()),
            'interior_cell_count':len(cells),'interior_slack_lower':str(min(g.arb(z[2]).lower() for z in cells)),
            'interior_degree_cells':cells,'side_cap_excess_per_mixing_weight_lower':str(end_sine.lower()),
            'middle_cap_excess_lower':str(mid_margin.lower()),'permitted_dihedral_poses':admitted,'excluded_dihedral_poses':excluded}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.write_text(json.dumps(verify(),indent=2)+'\n')
