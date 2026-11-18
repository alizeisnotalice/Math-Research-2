"""Certify the regular pentagon placement transition and stability bounds."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sys
import json
from flint import arb, ctx, __version__

ctx.prec=192
def q(x):
    x=F(x)
    return arb(x.numerator)/x.denominator
def angle(x):
    return q(x)*arb.pi()/180
def co(x):
    return angle(x).cos()
def si(x):
    return angle(x).sin()
def bounds(x):
    def asrat(z):
        a,k=z.man_exp()
        return str(F(int(a))*F(2)**int(k))
    return [asrat(x.lower()),asrat(x.upper())]
def run():
    r=1/(2*co(18))
    py=2*(r*co(6)-q('.5'))/arb(3).sqrt()
    peak=r+py
    X=r*co(24)+py*si(210)
    Y=r*si(24)+py*co(210)
    threshold=2*(Y/X).atan()*180/arb.pi()-30
    assert 14<threshold<15
    def A_arc(x):
        return q('.5')-r*co(36-x)+py*co(x)
    def B_arc(x):
        return q('.5')-r*co(12+x)+py*si(210+x)
    def B_vertex(x):
        return X*co(x)+Y*si(x)-q('.5')
    def C_arc(x):
        return q('.5')-r*co(12-x)+py*si(330+x)
    # Arc-projection signs and concavity in the three sectors in the proof.
    assert peak*co(18)>q('.5') and -r*si(6)-py*si(30)<0
    assert B_arc(0)<0 and B_arc(6)<0 and B_vertex(14)<0<B_vertex(15)
    assert C_arc(0)<0
    S=angle(0).union(angle(30))
    assert r*(angle(36)-S).cos()-py*S.cos()>0
    assert r*(angle(12)+S).cos()-py*(angle(210)+S).sin()>0
    assert r*(angle(12)-S).cos()-py*(angle(330)+S).sin()>0
    assert X*S.cos()+Y*S.sin()>0
    # Positive and negative slant limits on the closed physical sigma range.
    S=angle('1.50').union(angle('1.54'))
    Af=peak*S.cos()-q('.5')
    Bf=q('.5')-r*(angle(12)+S).cos()+py*(angle(210)+S).sin()
    Cf=q('.5')-r*(angle(12)-S).cos()+py*(angle(330)+S).sin()
    slant_gap=min(Af.lower(),(-Bf).lower(),(-Cf).lower())
    eps=q('.0001')
    toler=F('0.34')
    root_gap=r*(2*co(12)-1)*co(6-toler)-q('.5')-3*eps
    assert root_gap>0
    change=r*si(18)*angle(toler)
    E=change+(2/arb(3).sqrt())*(eps+change)
    assert slant_gap>E+q('.0227')
    nominal=[peak*co(F('1.54')+k) for k in (0,10,15)]
    robust=[a-E for a in nominal]
    return dict(width_one_radius=bounds(r),fixed_translation_y=bounds(py),fixed_top_vertex_y=bounds(peak),
                transition_degrees=bounds(threshold),cap_exclusion_margin=bounds(slant_gap),
                original_six_independent_maximum_shift='1/10000',two_caps_independent_maximum_shift='227/10000',
                orientation_error_degrees=str(toler),root_safety=bounds(root_gap),hausdorff_error=bounds(E),
                fixed_top_support_center_ten_fifteen=list(map(bounds,nominal)),
                robust_top_support_center_ten_fifteen=list(map(bounds,robust)),precision=ctx.prec,
                python=sys.version,flint=__version__,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    j=run()
    Path(sys.argv[1]).write_text(json.dumps(j,indent=2)+'\n')
    for k in ('transition_degrees','fixed_top_vertex_y','cap_exclusion_margin','hausdorff_error'):
        print(k,float(F(j[k][0])))
    for k in ('fixed_top_support_center_ten_fifteen','robust_top_support_center_ten_fifteen'):
        print(k,[float(F(a[0])) for a in j[k]])
