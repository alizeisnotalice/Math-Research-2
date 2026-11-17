#!/usr/bin/env python3
"""Exact fixed RM(3,6) obstruction and complete complement-cover audit."""
import argparse
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys
import time


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();start=time.monotonic()
    linear=set()
    for u,v,w in combinations(range(1,64),3):
        points={0,u,v,w,u^v,u^w,v^w,u^v^w}
        if len(points)==8:linear.add(tuple(sorted(points)))
    flats={sum(1<<(base^x) for x in L) for L in linear for base in range(64)}
    assert len(linear)==1395 and len(flats)==11160
    monomials=[14,21,35]
    H=sum(1<<x for x in range(64) if sum((x&m)==m for m in monomials)%2)
    assert H.bit_count()==16 and all(m.bit_count()==3 for m in monomials)
    assert not any(F&H==F for F in flats)
    complement=((1<<64)-1)^H
    eligible=sorted(F for F in flats if F&complement==F)
    by_point={p:[F for F in eligible if F>>p&1] for p in range(64) if complement>>p&1}
    nodes=0
    @lru_cache(None)
    def cover(remaining):
        nonlocal nodes
        nodes+=1
        if not remaining:return ()
        options=min(([F for F in by_point[p] if F&remaining==F]
                     for p in by_point if remaining>>p&1),key=len)
        for F in options:
            answer=cover(remaining^F)
            if answer is not None:return (F,)+answer
        return None
    answer=cover(complement)
    assert answer is None
    result={'status':'passed','arithmetic':'integer/F2','tolerance':0,
            'ANF_monomial_masks':monomials,'support':[x for x in range(64) if H>>x&1],
            'degree':3,'weight':16,'RM_minimum_weight':8,
            'linear_three_spaces':len(linear),'affine_three_flats':len(flats),
            'minimum_flats_in_support':0,'affine_three_flats_in_complement':len(eligible),
            'complement_six_flat_cover_exists':False,'exact_cover_complete':True,
            'exact_cover_states':nodes,'elapsed_seconds':time.monotonic()-start,
            'scope':'Refutes degree/weight-only splitting; NOT a partial-packing completion counterexample',
            'python':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
