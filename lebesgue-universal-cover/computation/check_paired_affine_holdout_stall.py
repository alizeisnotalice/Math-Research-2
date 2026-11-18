"""Certify five-body area on the first thin unresolved box in the blind comparison."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, ctx

KEY=[1205,151326]
FINE=2**44


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(a):
    sys.path.insert(0,str(a.comp.resolve()))
    from certify_five_harmonic import OLD,placement_area,sqrt_upper,validate_geometry
    from certify_nine import ar
    data=json.loads(a.audit.read_bytes())
    assert F(data['target'])==F('.83744') and data['node_limit']==1500000 and data['depth_limit']==135
    row=next(r for r in data['results'] if r['key']==KEY)
    assert row['paired'][:2]==[0,55972] and row['paired'][4:6]==[135,2]
    lo,hi=row['first_unresolved_box']
    root_lo,root_hi=row['root']
    assert all(l<=x<y<=h for l,x,y,h in zip(root_lo,lo,hi,root_hi,strict=True))
    center=[F(x+y,2*FINE) for x,y in zip(lo,hi,strict=True)]
    half=[F(y-x,2*FINE) for x,y in zip(lo,hi,strict=True)]
    ctx.prec=160
    pi=arb.pi()
    validate_geometry(pi)
    assert pi<ar(F(22,7))
    displacement=[]
    for n,j,t,_,slope in OLD:
        assert 1/(2*(pi/(2*n)).cos())<ar(F(3,5))
        assert sum((abs(center[k])+half[k])**2 for k in (t,t+1))<F(2,5)**2
        rotation=slope*half[j] if j is not None else F(0)
        displacement.append(rotation+sqrt_upper(half[t]**2+half[t+1]**2))
    delta=max(displacement)
    area=placement_area(center,pi,k=2048)
    center_lower=F(area['lower'])
    box_lower=center_lower-F(44,7)*delta-F(22,7)*delta**2
    assert box_lower>F('.837551707765')>F(data['target'])
    return dict(key=KEY,audit_sha256=digest(a.audit),source_sha256=digest(Path(__file__)),
                position_denominator=FINE,samples_per_arc=2048,box=[lo,hi],
                center_lower=str(center_lower),delta=str(delta),box_lower=str(box_lower),
                box_lower_at_1e12=str(F(int(box_lower*10**12),10**12)))


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    for name in ('comp','audit'):
        parser.add_argument(name,type=Path)
    print(json.dumps(run(parser.parse_args()),separators=(',',':')))
