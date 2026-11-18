#!/usr/bin/env python3
"""Replay all numerical support branches and exhaust small abstract configurations."""

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path

from generate_support import LEFT,RIGHT
from support_cube import closure,saturated_pairs


def bits(n,m):
    return {i for i in range(n) if (m>>i)&1}


def abstract_test():
    checks=Counter()
    for a,b in product(range(1,4),repeat=2):
        sa0,sb0=set(range(a+1)),set(range(b+1))
        pa,pb={0,a},{0,b}
        if a==b:
            try:saturated_pairs(a,b,((a,0),(0,b)))
            except ValueError:
                checks["double_saturation"]+=1
                continue
            raise AssertionError("double saturation")
        ban0=saturated_pairs(a,b,((a,0),(0,b)))
        all_pairs=set(product(sa0,sb0))
        for hole in (None,*range(1,a+b)):
            for requested in (None,*range(1,a+b+1)):
                required=[] if requested is None else [
                    {(i,j) for i,j in all_pairs if i+j==requested and i>0 and j>0}]
                zeros=set() if hole is None else {hole}
                state=closure(a,b,pa,pb,sa0,sb0,{a,b},zeros,ban0,required)
                checks["abstract_cases"]+=1
                for u,v in product(range(1<<(a+1)),range(1<<(b+1))):
                    S,T=bits(a+1,u),bits(b+1,v)
                    pairs=set(product(S,T))
                    if not pa<=S or not pb<=T or pairs&ban0 or any(i+j==hole for i,j in pairs):
                        continue
                    if not all(pairs&test for test in required):continue
                    assert state is not None
                    l,r,L,R,one,zero=state
                    attained={i+j for i,j in pairs}
                    assert l<=S<=L and r<=T<=R and one<=attained and not zero&attained
                    checks["feasible_abstract_supports"]+=1
    assert checks["feasible_abstract_supports"]>0
    return checks


def current(data):
    groups=Counter()
    fields=("A_positive","B_positive","A_allowed","B_allowed","one","zero")
    for nrow in data["degrees"]:
        N=nrow["degree"]
        for row in nrow["directions"]:
            a,b=row["degrees"]
            ban=saturated_pairs(a,b,((a,0),(0,b)))
            for rec in row["branches"]:
                pa,pb=set(LEFT)|{a},set(RIGHT)|{b}
                sa,sb=set(range(8,a+1))|pa,set(range(8,b+1))|pb
                for k,side in enumerate(rec["shape"],8):
                    if side in "AX":pa.add(k)
                    if side in "BX":pb.add(k)
                    if side=="A":sb.discard(k)
                    if side=="B":sa.discard(k)
                found=closure(a,b,pa,pb,sa,sb,set(range(32))-{7}|{N},{7},ban,
                              [{(i,31-i) for i in range(1,31)}])
                assert found is not None
                assert tuple(sorted(x) for x in found)==tuple(rec[k] for k in fields)
                groups[N]+=1
    return groups


def main():
    if not __debug__:raise SystemExit("assertions required")
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    checks=abstract_test()
    covered=current(json.loads(a.input.read_bytes()))
    result=dict(status="verified",tiny_cases=checks,branches_rebuilt=covered)
    with a.out.open("x") as f:f.write(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":main()
