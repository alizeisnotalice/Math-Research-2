"""Replay the sharp endpoint slope for a second normalized hexagon orientation."""
import argparse
import json
from pathlib import Path
from flint import arb, ctx


def run():
    ctx.prec = 192
    L = arb.pi()/3
    theta = L/2
    tau = (L/4).tan()
    s = tau/(1+(2+tau*tau).sqrt())
    projection = (1+theta.cos())/2
    rise = s+theta.sin()/2
    gamma = arb.atan2(rise,projection)
    # The two formulas are algebraically equal; the interval checks also locate the switch.
    assert (projection**2+rise**2).contains(1)
    assert (s-((1-projection**2).sqrt()-theta.sin()/2)).contains(0)
    assert 0<gamma<theta and 2*tau < arb(3)/5
    assert (3*s>arb('.3295303066473057')) and (3*s<arb('.3295303066473059'))
    def minimum_ramp(z):
        d = (2*z**2-z**4).sqrt()
        return 2*((z+d)/(1-z**2)).atan()

    examples = []
    for raw in ('.012','.1','.2','.3','.329'):
        z=arb(raw)/3
        length=L-2*minimum_ramp(z)
        assert length>0
        examples.append(dict(slope=raw,longest_zero_interval_degrees=str(180*length/arb.pi())))
    return dict(precision=ctx.prec,critical_slope=str(3*s),normalized_tangent=str(s),
                interior_root_degrees=30,first_curvature_switch_degrees=str(180*gamma/arb.pi()),
                extremal_unit_arcs=9,maximum_interval_examples=examples)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    opt=parser.parse_args()
    result=run()
    opt.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['critical_slope'],result['first_curvature_switch_degrees'])
