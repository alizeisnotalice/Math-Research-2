#!/usr/bin/env python3
"""Independent direct-character verification of saved minimum witnesses."""
import argparse
import json
from pathlib import Path


def direct(v):
    return [sum(y*(-1 if (a & x).bit_count()%2 else 1)
                for x,y in enumerate(v)) for a in range(len(v))]


def span_size(s):
    span={0}
    for a in s:
        span |= {x ^ a for x in span}
    return len(span)


def check(w,d):
    v=w['truth']
    assert all(x in (-1,1) for x in v)
    a=direct(v)
    assert a==w['walsh_numerators']
    assert sum(x!=0 for x in a)==w['s']
    assert span_size([i for i,x in enumerate(a) if x])==1 << d
    b=direct(w['zero_one_truth'])
    assert all(b[1 << j] for j in range(d))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('source',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    data=json.loads(a.source.read_text())
    count=0
    for row in data['exhaustive']:
        assert sum(x['count'] for x in row['histogram'])==row['functions']
        for d,w in row['minimum_by_rank'].items():
            check(w,int(d));count+=1
    w=json.loads(json.dumps(data['exhaustive'][-1]['minimum_by_rank']['4']))
    w['walsh_numerators'][0]+=1
    try:
        check(w,4)
    except AssertionError:
        rejected=True
    else:
        raise AssertionError('corrupted coefficient accepted')
    out={'status':'passed','direct_character_witnesses':count,
         'tampered_coefficient_rejected':rejected,
         'scope':'minimum witnesses and ledger counts; not a second exhaustive search'}
    a.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))


if __name__=='__main__':
    main()
