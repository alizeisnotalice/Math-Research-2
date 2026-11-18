"""Transfer one certified EH polygonal gain and support allowance to all members of an AH box."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx


def format_area(lo, hi):
    places=39
    scale=10**places
    endpoints=((scale*lo).__floor__(),(scale*hi).__ceil__())
    return [f'{ticks//scale}.{ticks%scale:0{places}d}' for ticks in endpoints]


def run(ah_path, eh_path, output):
    ctx.prec=160
    a,b=ah_path.read_bytes(),eh_path.read_bytes()
    A,B=json.loads(a),json.loads(b)
    assert A['complete'] and B['AH_sha256']==sha256(a).hexdigest()
    assert arb(B['whole_rectangle_AH_clearance'])>0
    el,eu=map(F,B['EH_gain_over_R126_exact_interval'])
    strong_l,strong_u=map(F,A['strong_member_area_exact'])
    weak_l,weak_u=map(F,A['weak_member_area_exact'])
    assert 0<el<=eu and strong_l<=strong_u<weak_l<=weak_u
    assert tuple(map(F,B['area_exact_interval']))==(strong_l-eu,strong_u-el)
    u=arb(A['uniform_cut_line_margin_lower'])
    cost=arb(B['AH_local_chord_transfer_cost_upper'])
    assert u-cost>0
    endpoints=strong_l-eu,weak_u-el
    out=dict(complete=True,all_members_area_exact_interval=list(map(str,endpoints)),
             all_members_area_decimal_interval=format_area(*endpoints),
             strongest_area_decimal_interval=B['area_decimal_interval'],
             weak_member_area_decimal_interval=format_area(weak_l-eu,weak_u-el),
             uniform_post_EH_slack_lower=str((u-cost).lower()),
             AH_family_receipt_sha256=sha256(a).hexdigest(),
             EH_certificate_sha256=sha256(b).hexdigest())
    output.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    for name in ('ah','eh','output'):
        cli.add_argument(name,type=Path)
    args=cli.parse_args()
    run(args.ah,args.eh,args.output)
