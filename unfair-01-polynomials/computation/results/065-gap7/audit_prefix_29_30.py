"""Cross-check sparse prefix expansions against direct integral-ring convolution."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from random import Random
import json

from prefix_29_30 import Box,SCALE,T,canon,exact_system,model,slots
from verify_degree_19 import element


def digest(path):return sha256(path.read_bytes()).hexdigest()


def powers(row,coordinates):
    y=element(0)
    for monomial,v in row.items():
        weight=1
        for index in monomial:
            weight*=coordinates[index]
        y+=weight*v
    return y


def main():
    if not __debug__:raise SystemExit('requires assertions')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--prefix',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    options=p.parse_args()
    cases=json.loads(options.prefix.read_bytes())
    count=Counter()
    for degree in cases:
        for idx,r in enumerate(degree['late']):
            if idx%19:continue
            j=degree['degree']
            shape=canon(r['shape'])
            unknown=slots(j,shape)
            if not unknown:continue
            n=len(unknown)
            equations=exact_system(j,shape)
            for k in range(4):
                numeric=[((idx+7*k+3*i+i*k)%13)-6 for i in range(n)]
                raw=[(r,powers(f,numeric)) for r,f in equations]
                assert raw==model(T,j,shape,numeric)[0]
                count[str(j)+'_assignments']+=1
            count[str(j)+'_shapes']+=1
            count[str(j)+'_quadratic']+=int(any(any(len(m)==2 for m in f) for _,f in equations))
    rng=Random(29030)
    for k in range(4096):
        x,y=sorted(rng.randrange(-2*SCALE,2*SCALE+1) for _ in range(2))
        a,b=sorted(rng.randrange(-2*SCALE,2*SCALE+1) for _ in range(2))
        result=Box(x,y)*Box(a,b)
        values=[Fraction(h*i,SCALE**2) for h in (x,y) for i in (a,b)]
        assert Fraction(result.lo,SCALE)<=min(values)<=max(values)<=Fraction(result.hi,SCALE)
        count['rational_interval_product_cases']+=1
    data={'status':'verified','counts':dict(count),'seed':29030,
          'input_sha256':digest(options.prefix),'prefix_source_sha256':digest(Path(__file__).with_name('prefix_29_30.py'))}
    with options.out.open('x') as f:f.write(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(json.dumps(data,sort_keys=True))


if __name__=='__main__':main()
