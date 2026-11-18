"""Exact dihedral support signs, and curvature/slack budgets for the two unequal caps."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from flint import arb,ctx

def run():
    ctx.prec=192
    def six(e):
        return [e[j//2] if j%2==0 else -e[((j+3)%6)//2] for j in range(6)]
    cases=[]
    for hand,targets,shift in [(1,([2,2,-1],[5,-1,-1]),2),(-1,([2,2,-1],[-1,5,-1]),-2)]:
        a,b=map(six,targets)
        rows=[[j,(j+shift)%6,a[j],b[(j+shift)%6]] for j in range(6)]
        assert all(x*y<0 for j,k,x,y in rows)
        cases.append(dict(hand=hand,rows=rows))
    # Positive selector when s3=s5: use the residues of the low support signs.
    for mask in range(8):
        clear={j if (mask>>j)&1 else j+3 for j in range(3)}
        assert any(j in clear and (j+2)%6 in clear for j in range(6))
    w=F(1,40);eps=F(1,200000)
    # t variable on [-1,1]. The profile is constant on [-1/4,1/4], and varies
    # by quintic smoothstep outside; |b'| <= 5/2 and |b''| <= 12.
    radloss=eps*(40+66/w+53/w**2)
    assert radloss<F(1,2)
    a=lambda q:arb(q.numerator)/q.denominator
    rad=arb.pi()/180
    smin,smax=F('1.50'),F('1.54')
    assert 2*a(w)<2*a(smin)*rad
    assert a(smax-smin)*rad/2<a(w)/4
    assert 3*a(smax)*rad<arb.pi()/2
    slope=a(eps)*(3*a(smax)*rad).cos()/(2*a(w))
    assert slope>arb('.00009967')
    d=F('.04')
    rep=slope*a(d)*rad
    assert rep>arb('.00000006958')
    return dict(complete=True,scope=['1.50','1.54'],sign_rows=cases,phase_centers=['pi/6+tau','pi/6-tau'],
        profile_width_radians=str(w),epsilon=str(eps),curvature_radius_lower=str(F(1,2)-radloss),
        support_violation_per_radian_lower=str(slope.lower()),cap_gap_004_degrees_lower=str(rep.lower()),
        selector_orientation_checks=8,verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args()
    r=run();a.output.write_text(json.dumps(r,indent=1)+'\n');print(r['curvature_radius_lower'],r['cap_gap_004_degrees_lower'])
