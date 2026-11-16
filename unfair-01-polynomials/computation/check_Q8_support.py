"""Exact support checks used in the Q8 L>=9 proof."""
S = {0, 1, 2, 4, 5, 8}
T = {0, 1, 2, 4, 5}

for L in range(9, 201):
    sums = [i + L*j for i in S for j in S]
    assert len(sums) == len(set(sums)) == 36
    assert all({(s+h) % L for s in S} != S for h in range(1, L))

assert all({(t+h) % 9 for t in T} - S for h in range(1, 9))
print({"status": "passed", "L_checked": [9, 200],
       "Q8_support_size": 36, "Q8_nonzero_stabilizers": 0,
       "L9_T_to_S_nonzero_inclusions": 0})
