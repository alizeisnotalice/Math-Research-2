"""Certify two AH chords meeting on the retained straight region."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import ctx
from wider_chart import geometry, green, bracket
from certify_ah_margin_envelope import coefficients
from certify_ah_three_line_cut import dyadic, fit_archived_witnesses, format_area


def prepare(archive, lines):
    ctx.prec = 128
    extra = [tuple(map(F, line)) for line in lines]
    (p1,d1),(p2,d2) = extra
    assert F('.419')<p1<F('.430') and 0<p2<p1 and 0<d1<d2<F('.004')
    ar, geo = geometry(archive, extra)
    return ar, geo, extra


def area_and_witnesses(archive, ar, extra):
    import audit_optimized_eh_chord as g
    from audit_eh_continuous_placement import polygon
    from certify_eh_polyline_obstruction import clip_geometry

    ctx.prec = 160
    (p1,d1),(p2,d2) = extra
    A=ar['A']
    old=coefficients(ar)[:3]
    slope1,slope2 = (A(d).tan() for p,d in extra)
    intercept1,intercept2 = A(p1)*slope1,A(p2)*slope2
    assert intercept2>intercept1 and slope2>slope1
    split=(intercept2-intercept1)/(slope2-slope1)
    assert 0<split<A(p2)<A(p1)
    h=intercept1-slope1*split
    assert h>0
    # The longer, shallower chord already dominates the three anchored chords.
    assert all(intercept1-slope1*w>a-b*w for b,a,_ in old for w in (-A(ar['pbound']),A(p1)))
    tail = ar['length']**2*h*A(p1-p2)/2
    assert tail>0

    g.SIGMA=A('1.5497')*g.PI/180
    inner,_=polygon('0.0005522',32)
    dirs,offsets,*_=clip_geometry(['0.000555','0.0005522','0.0075'])
    inner.extend((f'R126_EH_{i}',u,c) for i,(u,c) in enumerate(zip(dirs,offsets)))
    # The triangle changing the tail has vertices at the kink and the two pivots.
    z,D,n,L = (ar[k] for k in ('z','D','normal','length'))
    kink=g.plus(g.plus(z,g.scale(split,D)),g.scale(L*h,n))
    pivots=[g.plus(z,g.scale(A(p),D)) for p in (p2,p1)]
    clearance=min((c-g.dot(u,x)).lower() for name,u,c in inner if name!='Ta-Z' for x in [kink]+pivots)
    assert clearance>0
    witnesses=fit_archived_witnesses(archive,ar,inner)
    cap_lo,cap_hi,shape=green(archive,p2,d2)
    lo_tail,hi_tail=map(dyadic,(tail.lower(),tail.upper()))
    shape.update(casting_kink_w=str(split),positive_height_at_kink=str(h),
                 added_straight_tail_triangle_exact=[str(lo_tail),str(hi_tail)],
                 old_boundary_tail_clearance_lower=str(clearance))

    def read(name):return json.loads((archive/name).read_bytes())
    gibbs=read('gibbs-1p5497-refined-18-128.json')
    eh=read('eh-two-line-cover-certificate.json')['summary']
    assert F(eh['sigma_degrees'])==F(gibbs['sigma_degrees'])==F('1.5497')
    assert F(eh['delta'])==F('.000555') and eh['AH_relocation_after_root'] and eh['delta_is_enclosing_placement_depth']
    a,b=map(F,gibbs['area_interval'])
    el,ehigh=map(F,eh['cut_gain_area_interval_exact'])
    lo,hi=a-ehigh-cap_hi-hi_tail,b-el-cap_lo-lo_tail
    assert lo<hi
    return lo,hi,shape,witnesses


def run(archive, lines, out):
    ar,geo,extra=prepare(archive,lines)
    lo,hi,shape,witnesses=area_and_witnesses(archive,ar,extra)
    ctx.prec=128
    stats, verifier = bracket(ar,None)
    result={'lines':lines,'area_exact_interval':[str(lo),str(hi)],
            'area_decimal_interval':format_area(lo,hi),'shape':shape,
            'R150_R156_witnesses_at_identity':witnesses,'bins':stats,
            'geometry_sha256':sha256(geo.encode()).hexdigest(),
            'fine_verifier_sha256':sha256(verifier.encode()).hexdigest(),
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':sys.version,'flint':ar['flint_version']}
    out.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    for name in ('archive','input','output'):
        cli.add_argument(name,type=Path)
    opt=cli.parse_args()
    sys.path.insert(0,str(opt.archive.resolve()))
    run(opt.archive.resolve(),json.loads(opt.input.read_bytes()),opt.output)
