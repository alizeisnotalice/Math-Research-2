#!/usr/bin/env python3
"""Exact index audit for R(x^a)R(x^9), gcd(a,9)=1, 2<=a<=8."""
from __future__ import annotations

import json
from math import gcd
from pathlib import Path

S = (0, 1, 2, 3, 4, 5, 8)
AS = tuple(a for a in range(2, 9) if gcd(a, 9) == 1)

def pairs(x: set[int], y: set[int], total: int) -> list[tuple[int, int]]:
    return [(i, total-i) for i in sorted(x) if total-i in y]

def main() -> None:
    rows = []
    for a in AS:
        target = [a*i + 9*j for i in S for j in S]
        assert len(set(target)) == 49
        carry = []
        for i in range(17):
            for j in range(17):
                same = [(r,s) for r in S for s in S if a*r+9*s == a*i+9*j]
                if same and (i,j) != same[0]: carry.append(((i,j),same[0]))
        expected = [((i+9,j-a),(i,j)) for i in range(6) for j in S if j >= a]
        assert carry == expected
        ii0x, ii0y = {0,1,2,5}, {0,3,8}
        ii1x, ii1y = {0,1,2,5,8}, {0,3}
        assert not pairs(ii0x,ii0y,11)
        assert all(not pairs(ii1x,ii1y,q) for q in (9,10,12,13,14))
        rows.append({"a":a,"carry_count":len(carry),"target_injective":True,
                     "high_v_degree_upper_bound":8-a,"a_in_S":a in S,
                     "II_required_zero_rows_verified":True})
    out={"status":"passed_exact_index_audit","a_values":list(AS),"S":list(S),"records":rows,
         "boundary":"Finite support ledger only; the proof is 024-R8-all-binary-substitutions.md."}
    print(json.dumps(out,indent=2))

if __name__ == '__main__': main()
