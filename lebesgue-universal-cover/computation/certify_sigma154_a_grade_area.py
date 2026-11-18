"""Certify the true area of the .429/.009 AH geometry with either proved EH profile."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_three_sigma_area_budget import fixed_ah,setup,pack,subtract,ah_clearance,exact,U0


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(a):
    comp=a.comp.resolve()
    sys.path.insert(0,str(comp))
    import audit_optimized_eh_chord as g
    ctx.prec=160
    first=json.loads(a.single.read_bytes())
    chain=json.loads(a.chain.read_bytes())
    assert chain['single_sha256']==digest(a.single)
    assert first['sigma']==chain['sigma']=='1.54'
    base=tuple(map(F,first['base_area']['exact']))
    assert base==tuple(map(F,chain['original_baseline']['exact']))
    theta,pivot,far=F('.009'),F('.429'),F('.273')
    cap,w,slope=fixed_ah(comp,'1.54',pivot,theta,far)
    assert w[1]<far
    geom=setup('1.54',g)
    U,C3,T=(geom[k] for k in ('U','C3','T'))
    r=g.minus(U,C3)
    result=[]
    for label,d,data in (('single',F(first['EH_depth']),first),('eleven',F(chain['EH_deepest']),chain)):
        c,s=g.unit(exact(d))
        V=g.plus(C3,g.minus(g.scale(c,r),g.scale(s,g.rot(r))))
        clearance=ah_clearance(geom,V,T,g,pivot,theta)
        gain=tuple(map(F,data['EH_cut']['exact']))
        total=subtract(subtract(base,cap),gain)
        saving=subtract((U0-F('1e-6'),)*2,total)
        assert saving[0]>0
        row=dict(label=label,EH_outer_depth=str(d),EH_cut=pack(gain),AH_EH_clearance=pack(clearance),
                 resulting_area=pack(total),saving_below_A=pack(saving))
        if label=='eleven':
            row['strict_AH_margin_needed']=chain['transfer_cost']
        result.append(row)
    return dict(sigma='1.54',pivot=str(pivot),turn=str(theta),AH_far=str(far),
                AH_curve_crossing=pack(w),AH_crossing_derivative_lower=str(slope),
                base=pack(base),AH_cut=pack(cap),cases=result,
                single_sha256=digest(a.single),eleven_sha256=digest(a.chain),
                source_sha256=digest(Path(__file__)),
                proof_scope="True geometry and areas. New AH placement at .009 is supplied separately.")


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    for key in ('comp','single','chain'):
        parser.add_argument(key,type=Path)
    print(json.dumps(run(parser.parse_args()),separators=(',',':')))
