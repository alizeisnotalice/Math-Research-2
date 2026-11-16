"""Exact macro-coordinate audit for H_d(x^3) H_d(x^5).

Coordinates (r,e) represent u**e*v**r after u=x, v=x**d.
The script prints target binary collisions and every actual-to-canonical
alias for d=6..45.  It uses only finite integer sets, so it is an exact
accounting certificate rather than a positivity proof.
"""
from collections import Counter, defaultdict
from pathlib import Path
import json
import sympy as sp


def p_coords(a):
    return Counter({
        (0, 0): 1, (0, a): 1, (0, 2 * a): 1,
        (a, 0): 1, (a, -a): 1,
    })


def mul(left, right):
    ans = Counter()
    for (r, e), c in left.items():
        for (s, f), q in right.items():
            ans[r + s, e + f] += c * q
    return ans


TARGET = mul(p_coords(3), p_coords(5))
T = set(TARGET)
ACTUAL = mul(TARGET, TARGET)


def image(coords, d):
    ans = defaultdict(list)
    for r, e in coords:
        ans[r * d + e].append((r, e))
    return ans


def aliases(d):
    canonical = image(T, d)
    injective = all(len(items) == 1 for items in canonical.values())
    if not injective:
        return injective, []
    result = []
    for (r, e), coefficient in sorted(ACTUAL.items()):
        for s, f in canonical.get(r * d + e, []):
            if (r, e) != (s, f):
                result.append({
                    "actual": [r, e],
                    "canonical": [s, f],
                    "full_target_squared_multiplicity": coefficient,
                    "canonical_coefficient": TARGET[(s, f)],
                })
    return injective, result


def possible_parameters(left, right):
    ans = set()
    for r, e in left:
        for s, f in right:
            if r != s and (f - e) % (r - s) == 0:
                d = (f - e) // (r - s)
                if d > 0:
                    ans.add(d)
    return ans


def main():
    assert sum(TARGET.values()) == 25
    assert all(c == 1 for c in TARGET.values())
    target_rows = {r: sorted(e for rr, e in T if rr == r) for r in sorted({r for r, _ in T})}
    actual_rows = {r: sorted(e for rr, e in ACTUAL if rr == r) for r in sorted({r for r, _ in ACTUAL})}
    print("T rows:", target_rows)
    print("T+T rows:", actual_rows)
    target_collision_parameters = possible_parameters(T, T)
    alias_parameters = possible_parameters(set(ACTUAL), T)
    assert {d for d in target_collision_parameters if d >= 6} == {6}
    assert max(alias_parameters) == 14
    print("all positive target collision parameters:", sorted(target_collision_parameters))
    print("all positive actual-to-canonical alias parameters:", sorted(alias_parameters))
    data = {}
    for d in range(6, 46):
        injective, found = aliases(d)
        data[str(d)] = {"target_macro_injective": injective, "aliases": found}
        if not injective:
            collisions = {n: v for n, v in image(T, d).items() if len(v) > 1}
            print(f"d={d}: TARGET NONBINARY {collisions}")
        elif found:
            print(f"d={d}: " + "; ".join(
                f"{z['actual']}->{z['canonical']}" for z in found))
        else:
            print(f"d={d}: no alias")

    expected_tail_pairs = {
        11: {(3, 5), (6, 8)},
        12: {(6, 5)},
        13: {(3, 5)},
        14: {(3, 5), (6, 8)},
    }
    for d, expected in expected_tail_pairs.items():
        injective, found = aliases(d)
        assert injective
        assert {
            (entry["actual"][0], entry["canonical"][0]) for entry in found
        } == expected, (d, found)
    for d in range(15, 46):
        assert aliases(d) == (True, []), d
    expected_small_counts = {
        8: {(0, 3): 4, (3, 5): 3, (5, 8): 2},
        9: {(0, 3): 4, (3, 5): 4, (5, 8): 2, (6, 5): 2, (6, 8): 2},
        10: {(0, 3): 2, (3, 5): 6, (5, 8): 1, (6, 8): 4},
    }
    for d, expected in expected_small_counts.items():
        injective, found = aliases(d)
        assert injective
        counts = Counter((entry["actual"][0], entry["canonical"][0]) for entry in found)
        assert dict(counts) == expected, (d, counts)

    # Exact gcd audit: gcd(L3*L5,H_d(u^3)H_d(u^5)) is Phi_15 exactly
    # for d == 4 mod 5.  The second gcd used in the non-pure branch is
    # always one.
    u = sp.symbols("u")
    L3u = 1 + u**3 + u**6
    L5u = 1 + u**5 + u**10
    phi15 = sp.cyclotomic_poly(15, u)
    for residue in range(45):
        d0 = residue or 45
        h = lambda z0: 1 + z0 + z0**2 + z0**(d0 - 1) + z0**d0
        g = sp.factor(sp.gcd(L3u * L5u, sp.expand(h(u**3) * h(u**5))))
        assert g == (phi15 if d0 % 5 == 4 else 1), (d0, g)
        assert sp.gcd(L5u, h(u**5)) == 1, d0

    # Symbolic row identity for the proposed 3->5, 6->8 transport.
    L3, L5, U3, U5, D, E, z = sp.symbols("L3 L5 U3 U5 D E z")
    V = L3 * L5
    R = U3 * L5 + z * D
    S = L3 * U5 - D
    T6 = z * E
    W = U3 * U5 - E
    assert sp.expand(
        R * S * W - S**2 * T6 - W**2 * V
        - (L3 * E - U3 * D) * (L5 * W - z * U5 * S)
    ) == 0

    # d=12 has the sole reverse alias actual (6,-6) -> canonical (5,6).
    # Its coefficient c has the displayed sum-of-nonnegative-terms identity.
    c = sp.symbols("c")
    V12 = L3 * L5
    R12 = U3 * L5
    S12 = L3 * U5 - c * z
    T12 = c / z
    W12 = U3 * U5
    assert sp.expand(
        R12 * S12 * W12 - S12**2 * T12 - W12**2 * V12
        + c / z * (S12**2 + L5 * U3**2 * U5 * z**2)
    ) == 0
    output = Path(__file__).with_name("hd35_macro_carries_d6_d45.json")
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print("row syzygy passed exactly; wrote", output)
    print("tail assertions: d>=15 no carry; d=11,13,14 expected transport only;")
    print("pure-row gcd is Phi_15 exactly for d == 4 mod 5.")
    print("d=12 reverse-carry sum-of-nonnegative-terms identity passed exactly.")
    print("d=8/9/10 exact alias-type counts passed.")


if __name__ == "__main__":
    main()
