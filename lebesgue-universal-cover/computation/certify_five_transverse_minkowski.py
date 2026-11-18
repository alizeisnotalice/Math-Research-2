"""Certify transverse growth using the Minkowski sum of six independent contact-chain choices."""
import argparse
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
from flint import arb,ctx

def main():
    if not __debug__:raise SystemExit("Assertions required")
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("contact",type=Path)
    p.add_argument("--weights",type=Path,default=Path(__file__).with_name("five-transverse-orthant-weights.json.gz"))
    a=p.parse_args()
    contact=a.contact.resolve()
    trans=contact/"certify_five_body_contact_transverse.py"
    stat=contact/"certify_five_body_contact_stationary.py"
    sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    assert sha(trans)=="6460ce26be5f1f4ae5ec74c1aeda7adb4032958810f2c10ce2316b3ffac10956"
    assert sha(stat)=="c46a06d878962083f005692a66b7a546779e4ee0a1cb4d42c075e8b1a3ed3a4f"
    sys.path.insert(0,str(contact))
    from certify_five_body_contact_transverse import gradients,BOX,interval
    ctx.prec=192
    base,groups,_=gradients(*(interval(*x) for x in BOX))
    widths=list(map(len,groups));assert widths==[8,8,8,8,4,2]
    weights=json.loads(gzip.decompress(a.weights.read_bytes()))
    assert len(weights)==512
    minimum=arb(1)
    for sign,row in zip(product((-1,1),repeat=9),weights,strict=True):
        g=base.copy()
        assert list(map(len,row))==widths
        for group,z in zip(groups,row):
            assert all(type(x)is int and x>=0 for x in z) and sum(z)==2**34
            for j in range(9):
                g[j]+=sum((w*val[j] for w,val in zip(z,group)),arb(0))/2**34
        lower=min((x*k).lower()for x,k in zip(g,sign))
        assert lower>arb(42804)/10**8
        minimum=min(minimum,lower)
    print(json.dumps(dict(events=6,choices_per_event=widths,material_chains=32768,orthants=len(weights),
        exact_integer_weight_denominator=2**34,transverse_growth_floor="0.00042804",
        minimum_verified_component_lower=str(minimum),weights_sha256=sha(a.weights),script_sha256=sha(Path(__file__)),
        transverse_geometry_sha256=sha(trans),stationary_geometry_sha256=sha(stat)),indent=2))

if __name__=="__main__":main()
