"""Independently verify selected targets and their exact model prerequisites."""

import argparse

from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

from sympy import Poly, ZZ

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "062-residue-sieve"))

from prepare_gap7_degree56 import X
from prepare_two_residue_targets import decode, value
from verify_gap7_periods import certify


HERE=Path(__file__).resolve().parent


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()

    raw=(HERE/'inputs-gap7-degree56.json').read_bytes()
    inputs=json.loads(raw)
    source=(HERE/'ramified-masks-two-residues.txt').read_bytes()
    assert sha256(source).hexdigest()==inputs['source_mask_list_sha256']
    masks=list(map(int, source.splitlines()))
    full_set=set(masks)
    assert len(masks)==len(full_set)==425==inputs['matched_words']
    for mask in masks:
        c=decode(mask)
        assert value(c,-2,163)==value(c,-2,163,True)==0
        assert value(c,2,1279)==value(c,2,1279,True)==0
    jobs=[]
    expected=[]
    h_periods={int(p):v for p,v in inputs['h_reduction_periods'].items()}
    seen=set()
    for record in inputs['rejected_draws']:
        assert record['kind']=='period'
        mask=record['mask']
        assert mask in full_set and mask not in seen
        seen.add(mask)
        p=record['prime']
        degrees=record['degrees']
        period=lcm(*(d for d,_ in degrees))
        assert period % h_periods[p]
        jobs.append((tuple(reversed(decode(mask))),p,degrees))
        expected.append(period)
    for task in inputs['tasks']:
        mask=task['mask']
        assert mask in full_set and mask not in seen
        seen.add(mask)
        c=decode(mask)
        assert c==task['C']
        polynomial=Poly.from_list(c[::-1],X,domain=ZZ)
        assert polynomial.gcd(polynomial.diff()).degree()==0
        for p in (163,1279):
            pp=Poly(polynomial,modulus=p)
            assert pp.gcd(pp.diff()).degree()>0
        for p in h_periods:
            value_at_p=task['periods'][str(p)]
            assert value_at_p['period']%h_periods[p]==0
            jobs.append((tuple(reversed(c)),p,value_at_p['degrees']))
            expected.append(value_at_p['period'])
    assert len(seen)==inputs['attempted_draws']==255
    with ProcessPoolExecutor(max_workers=4) as pool:
        assert list(pool.map(certify,jobs))==expected
    report={'status':'verified','packet_sha256':sha256(raw).hexdigest(),
            'matched_words':len(masks),'period_tests':len(jobs),'targets':len(inputs['tasks'])}
    with args.out.open('x') as output:
        json.dump(report,output,indent=2)
        output.write('\n')
    print(json.dumps(report))


if __name__=='__main__':
    main()
