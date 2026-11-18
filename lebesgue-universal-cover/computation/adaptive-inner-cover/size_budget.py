"""Evaluate the Hausdorff and area transfers for an explicit point-count budget."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, ctx
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_inner_approximation import one_row, rational


def first_true(predicate):
    high = 4
    while not predicate(high):
        high *= 2
    low = 2
    while high-low > 1:
        mid = (low+high)//2
        if predicate(mid):
            high = mid
        else:
            low = mid
    return high


def compute(comp, output):
    ctx.prec = 256
    pi = arb.pi()
    cache = {}

    def budget(n):
        if n not in cache:
            cache[n] = arb(one_row(n)['finite_upper']), (pi/n).cos()/2
        return cache[n]

    def strict_above(x,y):
        if x > y:
            return True
        assert x < y, (x,y)
        return False

    def ceiling_for(x,n):
        u,r = budget(n)
        return (1+2*u/r)*x+pi*u**2

    def lower_for(x,n):
        u,r = budget(n)
        return (x-pi*u**2)/(1+2*u/r)

    orig = comp/'finite-b-core'/'host-area.json'
    tight = comp/'single-repair'/'single-repair.json'
    face = comp/'released-hexagon'/'difference.json'
    j = lambda p: json.loads(p.read_bytes())
    target = rational('.844082046521218581802948589770767704354')
    specs = {'round215':j(orig)['host_area_upper'],'round219':j(tight)['base_area'][1],
             'round223':j(face)['single']['upper']}
    upper_profiles = {}
    for name,value in specs.items():
        a = rational(value)
        n = first_true(lambda k: strict_above(target,ceiling_for(a,k)))
        assert ceiling_for(a,n) < target < ceiling_for(a,n-1)
        upper_profiles[name] = {'area':value,'first_sufficient_N':n,'at_N':str(ceiling_for(a,n)),
                                'at_previous':str(ceiling_for(a,n-1))}
    lower_source = comp/'six-input-083745'/'points'/'normal-gap-083745-audit.json'
    lower_data = j(lower_source)
    lower_fraction = Fraction(3352479567681343194888546017449042074380160550693261,
                              4003192723689732291403922205052321660928000000000000)
    assert str(lower_fraction) in json.dumps(lower_data)
    n = first_true(lambda k: strict_above(lower_for(rational(lower_fraction),k), rational('.83745')))
    assert lower_for(rational(lower_fraction),n)>rational('.83745')>lower_for(rational(lower_fraction),n-1)
    area4 = pi/8 + arb(3).sqrt()/4
    assert arb('.8257') < area4 < arb('.8258')
    rows = {'B_area_budget':upper_profiles,'four_point_candidate_area':str(area4),
            'lower_083745_transfer':{'general_source':str(lower_fraction),'first_sufficient_N':n,
                    'at_N':str(lower_for(rational(lower_fraction),n)),
                    'at_previous':str(lower_for(rational(lower_fraction),n-1))},
            'precision_bits':ctx.prec,'source_sha256':{str(p.relative_to(comp)):sha256(p.read_bytes()).hexdigest()
                                           for p in (orig,tight,face,lower_source)}}
    output.write_text(json.dumps(rows,indent=2)+'\n')
    print({k:v['first_sufficient_N'] for k,v in upper_profiles.items()},'lower',n)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('comp',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    assert __debug__
    compute(args.comp,args.output)
