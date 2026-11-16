#!/usr/bin/env python3
"""Exact support/row audit for R(x^a)R(x^10), a=3,7."""
from __future__ import annotations

import json
from pathlib import Path

S = (0, 1, 2, 3, 4, 5, 8)

def pairs(x: set[int], y: set[int], n: int) -> list[tuple[int,int]]:
    return [(i,n-i) for i in sorted(x) if n-i in y]

def main() -> None:
    records=[]
    for a in (3,7):
        target=[a*i+10*j for i in S for j in S]
        assert len(set(target)) == 49
        actual=[]
        for i in range(17):
            for j in range(17):
                hit=[(r,s) for r in S for s in S if a*r+10*s==a*i+10*j]
                if hit and hit[0] != (i,j): actual.append(((i,j),hit[0]))
        expected=[((i+10,j-a),(i,j)) for i in range(6) for j in S if j>=a]
        assert actual == expected
        # II after the actual zero K9 has removed Y8 and Y1, and after u6/u7.
        ii_x,ii_y={0,1,2,5,8},{0,2,3}
        assert pairs(ii_x,ii_y,9) == []
        assert all(pairs(ii_x,ii_y,n)==[] for n in (12,13,14,15))
        # If X5!=0, K7=0 is the legal reason that Y2 disappears.
        assert pairs(ii_x,ii_y,7) == [(5,2)]
        records.append({"a":a,"target_grid_injective":True,
                        "carry_count":len(actual),"high_v_degree_bound":8-a,
                        "II_K7_pairs":pairs(ii_x,ii_y,7),
                        "II_high_voids":True})
    out={"status":"passed_exact_support_audit","S":list(S),"records":records,
         "boundary":"Finite index audit only; 024-R8-all-binary-substitutions.md gives the all-coefficient proof."}
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
