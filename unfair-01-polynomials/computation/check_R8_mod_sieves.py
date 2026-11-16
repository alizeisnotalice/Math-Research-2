#!/usr/bin/env python3
"""Exact residue-support audit for the proposed R8 L=10/11 arguments."""
import json

E = {0,1,2,3,4,5,8}

def candidates(L, fixed_other):
    return sorted(r for r in E if all((r+s) % L in E for s in fixed_other))

out = {}
# II: A has low residues 0,1,2,8; B has 0,3.
for L in (10,11):
    out[f"II_L{L}"] = {
        "A_candidates_from_B_0_3": candidates(L, {0,3}),
        "B_candidates_from_A_0_1_2_8": candidates(L, {0,1,2,8}),
    }
# III, L=10: A has 0,1,8; B has 0,2,4.
out["III_L10"] = {
    "A_candidates_from_B_0_2_4": candidates(10, {0,2,4}),
    "B_candidates_from_A_0_1_8": candidates(10, {0,1,8}),
}
# I: with all E present in A, only B residue 0 survives for these L.
for L in (10,11,12):
    out[f"I_L{L}"] = {"B_candidates_from_full_A_E": candidates(L, E)}
print(json.dumps(out, ensure_ascii=False, indent=2))
