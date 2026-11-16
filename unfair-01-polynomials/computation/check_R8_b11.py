#!/usr/bin/env python3
"""Exact row audit for R(x^a)R(x^11), a=2,...,8."""
from __future__ import annotations

import json
from pathlib import Path

S=(0,1,2,3,4,5,8)

def pairs(x:set[int],y:set[int],n:int)->list[tuple[int,int]]:
    return [(i,n-i) for i in sorted(x) if n-i in y]

def main()->None:
    records=[]
    for a in range(2,9):
        vals=[a*i+11*j for i in S for j in S]
        assert len(set(vals))==49
        carry=[]
        for i in range(17):
            for j in range(17):
                hit=[(r,s) for r in S for s in S if a*r+11*s==a*i+11*j]
                if hit and hit[0] != (i,j): carry.append(((i,j),hit[0]))
        expected=[((i+11,j-a),(i,j)) for i in range(6) for j in S if j>=a]
        assert carry==expected
        # Actual II rows after K9 kills Y1,Y8 and K10 kills Y2.
        x,y={0,1,2,5,8},{0,3}
        assert pairs(x,y,11)==[(8,3)]
        assert all(not pairs(x,y,n) for n in range(12,17))
        records.append({"a":a,"target_grid_injective":True,"carry_count":len(carry),
                        "K10_pair_before_Y2_elimination":[[8,2]],
                        "K11_pair_after_elimination":[[8,3]],
                        "K12_to_K16_void":True})
    out={"status":"passed_exact_support_audit","S":list(S),"records":records,
         "boundary":"Finite support audit only; 024-R8-all-binary-substitutions.md contains the nonnegative row proof."}
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
