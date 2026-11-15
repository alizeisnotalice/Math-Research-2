#!/usr/bin/env python3
"""Pure graph check for the embedded S8 and 20-hole extension."""
from fractions import Fraction
from collections import deque
import json

S8 = {(2, 0), (0, 1), (1, 1), (2, 1)}
S8 |= {(x, 0) for x in range(3, 8)}
S8 |= {(0, y) for y in range(2, 8)}
P8 = Fraction(
    198947127410665164977798511882264205133,
    686832275881468514216635111007873246808,
)

def adj(v, n):
    x, y = v
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        w = x + dx, y + dy
        if 0 <= w[0] < n and 0 <= w[1] < n:
            yield w

def connected(S, n):
    if not S:
        return False
    seen = {next(iter(S))}
    q = deque(seen)
    while q:
        v = q.popleft()
        for w in adj(v, n):
            if w in S and w not in seen:
                seen.add(w)
                q.append(w)
    return seen == S

def components(R, n):
    ans = []
    left = set(R)
    while left:
        seed = left.pop()
        c = {seed}
        q = deque([seed])
        while q:
            v = q.popleft()
            for w in adj(v, n):
                if w in left:
                    left.remove(w)
                    c.add(w)
                    q.append(w)
        ans.append(c)
    return ans

def boundary(S, n):
    return {w for v in S for w in adj(v, n) if w not in S}

def two_step_pairs(S, n):
    bd = boundary(S, n)
    closed = S | bd
    return {(x, y) for x in bd for y in adj(x, n) if y not in closed}

def prefix(n, holes=False):
    V = {(x, y) for x in range(n) for y in range(n)}
    shift = n - 8
    Slocal = {(x + shift, y + shift) for x, y in S8}
    block = {(x + shift, y + shift) for x in range(8) for y in range(8)}
    S = (V - block) | Slocal
    H = {(2 * i + 1, 2 * j + 1) for i in range(5) for j in range(4)}
    if holes:
        if not H <= S:
            raise AssertionError("a claimed hole was not in the base prefix")
        S -= H
    else:
        H = set()
    return V, S, H

rows = []
for n in (18, 19, 32):
    V, S0, _ = prefix(n, holes=False)
    V, S, H = prefix(n, holes=True)
    R0, R = V - S0, V - S
    bd0, bd = boundary(S0, n), boundary(S, n)
    pairs0, pairs = two_step_pairs(S0, n), two_step_pairs(S, n)
    cs = components(R, n)
    row = {
        "n": n,
        "base": {
            "S_size": len(S0), "R_size": len(R0), "boundary": len(bd0),
            "connected_S": connected(S0, n), "two_step_pairs": len(pairs0),
        },
        "holes": {
            "count": len(H), "S_size": len(S), "R_size": len(R),
            "boundary": len(bd), "connected_S": connected(S, n),
            "two_step_pairs": len(pairs),
            "H_independent": all(w not in H for h in H for w in adj(h, n)),
            "H_no_edge_to_R8": all(w not in R0 for h in H for w in adj(h, n)),
            "R_component_sizes": sorted(map(len, cs)),
        },
    }
    assert row["base"] == {
        "S_size": n * n - 49, "R_size": 49, "boundary": 14,
        "connected_S": True, "two_step_pairs": 12,
    }
    assert row["holes"]["R_size"] == 69
    assert row["holes"]["boundary"] == 34
    assert row["holes"]["connected_S"]
    assert row["holes"]["two_step_pairs"] == 12
    assert row["holes"]["H_independent"]
    assert row["holes"]["H_no_edge_to_R8"]
    assert row["holes"]["R_component_sizes"] == [1] * 20 + [2, 47]
    rows.append(row)

pnew = Fraction(49 * 48, 69 * 68) * P8
assert 34 * pnew == Fraction(392, 23) * P8
out = {
    "status": "passed",
    "rows": rows,
    "p8": [P8.numerator, P8.denominator],
    "p_new": [pnew.numerator, pnew.denominator],
    "boundary_times_p_new": [(34 * pnew).numerator, (34 * pnew).denominator],
    "boundary_times_p_new_decimal": float(34 * pnew),
}
print(json.dumps(out, indent=2, sort_keys=True))
