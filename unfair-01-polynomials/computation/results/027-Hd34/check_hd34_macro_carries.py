"""Exact macro-coordinate carry audit for H_d(x^3) H_d(x^4).

Use u=x, v=x^d.  A coordinate (r,e) means u**e v**r.
The canonical target is P_3 P_4.  For every d in 11..35 this
enumerates every product coordinate in T+T whose one-variable image
equals that of a *different* canonical coordinate in T.  It also checks
whether the canonical map itself is injective.

All arithmetic is integer/set arithmetic.  This is an accounting
certificate, not a proof that nonnegative factors do or do not exist.
"""
from collections import Counter, defaultdict
import json
from pathlib import Path
import sympy as sp


def p_coords(a):
    # 1+u^a+u^(2a)+(1+u^(-a))v^a
    return Counter({(0, 0): 1, (0, a): 1, (0, 2 * a): 1,
                    (a, 0): 1, (a, -a): 1})


def mul_coords(left, right):
    ans = Counter()
    for (r, e), c in left.items():
        for (s, f), q in right.items():
            ans[r + s, e + f] += c * q
    return ans


TARGET = mul_coords(p_coords(3), p_coords(4))
T = set(TARGET)


def canonical_image(d):
    image = defaultdict(list)
    for coord in T:
        r, e = coord
        image[r * d + e].append(coord)
    return image


def aliases(d):
    image = canonical_image(d)
    injective = all(len(v) == 1 for v in image.values())
    if not injective:
        return injective, []
    actual = mul_coords(TARGET, TARGET)
    out = []
    for (r, e), multiplicity in sorted(actual.items()):
        n = r * d + e
        for canonical in image.get(n, []):
            if (r, e) != canonical:
                out.append({
                    "actual": [r, e],
                    "canonical": list(canonical),
                    "actual_multiplicity_if_full_target_squared": multiplicity,
                    "canonical_coefficient": TARGET[canonical],
                })
    return injective, out


def possible_positive_alias_parameters(left, right):
    """All positive integer d solving r*d+e=s*d+f for distinct layers."""
    possible = set()
    for r, e in left:
        for s, f in right:
            if r == s:
                continue
            numerator = f - e
            denominator = r - s
            if numerator % denominator == 0 and numerator // denominator > 0:
                possible.add(numerator // denominator)
    return possible


def human_rows(d, data):
    if not data:
        return "none"
    return "; ".join(
        f"({z['actual'][0]},{z['actual'][1]}) -> "
        f"({z['canonical'][0]},{z['canonical'][1]})"
        for z in data
    )


def main():
    assert sum(TARGET.values()) == 25
    assert all(c == 1 for c in TARGET.values()), TARGET
    target_rows = {r: sorted(e for rr, e in T if rr == r) for r in sorted({r for r, _ in T})}
    actual_support = set(mul_coords(TARGET, TARGET))
    actual_ranges = {
        r: [min(e for rr, e in actual_support if rr == r),
            max(e for rr, e in actual_support if rr == r)]
        for r in sorted({r for r, _ in actual_support})
    }
    print("canonical target rows:", target_rows)
    print("T+T exact row offset ranges:", actual_ranges)
    table = {}
    for d in range(6, 36):
        injective, data = aliases(d)
        table[str(d)] = {"target_macro_injective": injective, "aliases": data}
        if injective:
            print(f"d={d}: {human_rows(d, data)}")
        else:
            collisions = {
                n: coords for n, coords in canonical_image(d).items() if len(coords) > 1
            }
            print(f"d={d}: TARGET NOT BINARY; canonical collisions {collisions}")
    target_collision_d = possible_positive_alias_parameters(T, T)
    assert {d for d in target_collision_d if d >= 6} == {6, 8, 9, 12}, target_collision_d
    alias_d = possible_positive_alias_parameters(set(mul_coords(TARGET, TARGET)), T)
    assert max(alias_d) == 26, alias_d
    assert not {d for d in alias_d if d >= 27}
    allowed_pairs = {(3, 4), (6, 7)}
    surviving_parameters = defaultdict(set)
    for d in range(11, 27):
        if d == 12:
            continue
        _, data = aliases(d)
        pairs = {(item["actual"][0], item["canonical"][0]) for item in data}
        assert pairs <= allowed_pairs, (d, data)
        if d >= 11:
            for pair in pairs:
                surviving_parameters[pair].add(d)
    assert dict(surviving_parameters) == {
        (3, 4): {11, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 26},
        (6, 7): {11, 13, 14, 15, 16, 17, 19, 20, 23},
    }, surviving_parameters
    canonical_r7 = {e for r, e in T if r == 7}
    for d in range(17, 27):
        _, data = aliases(d)
        supplied_from_r6 = {
            item["canonical"][1]
            for item in data
            if item["actual"][0] == 6 and item["canonical"][0] == 7
        }
        assert supplied_from_r6 != canonical_r7, (d, supplied_from_r6)
    r6_offsets = {e for r, e in actual_support if r == 6}
    assert r6_offsets == {-6, -3, -2, 0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 16}
    assert {d for d in range(11, 27) if d in r6_offsets} == {12, 13, 16}
    # In the binary range d != 12, only d=13,16 can supply canonical (7,0)
    # from actual row 6, hence only they can have W=0.

    # Exact symbolic verification of the row syzygy used in the proof.
    L3, L4, U3, U4, D, E, z = sp.symbols("L3 L4 U3 U4 D E z")
    V = L3 * L4
    R = U3 * L4 + z * D
    S = L3 * U4 - D
    Trow = z * E
    W = U3 * U4 - E
    assert sp.expand(
        R * S * W - S**2 * Trow - W**2 * V
        - (L3 * E - U3 * D) * (L4 * W - z * U4 * S)
    ) == 0

    # The two Laurent/polynomial gcd statements are periodic in d modulo 36.
    u = sp.symbols("u")
    L3u = 1 + u**3 + u**6
    L4u = 1 + u**4 + u**8
    for residue in range(36):
        d0 = residue or 36
        h = lambda z0: 1 + z0 + z0**2 + z0**(d0 - 1) + z0**d0
        c = sp.expand(h(u**3) * h(u**4))
        assert sp.gcd(L3u * L4u, c) == 1, residue
        assert sp.gcd(L4u, h(u**4)) == 1, residue

    output = Path(__file__).with_name("hd34_macro_carries_d6_d35.json")
    output.write_text(json.dumps(table, indent=2, sort_keys=True) + "\n")
    print(f"wrote {output}")
    print("Canonical target collision parameters d>=6: {6, 8, 9, 12}.")
    print("All actual-to-canonical alias parameters have d<=26; hence no carry for d>=27.")
    print("Surviving d>=11 layer pairs:", dict(surviving_parameters))
    print("Actual r=6 offsets:", sorted(r6_offsets))
    print("In the binary d>=11 range, W=0 is possible only at d=13 or d=16.")
    print("For every d=17..26, the r=6 aliases do not cover the whole canonical r=7 row.")
    print("Row syzygy and period-36 gcd checks passed exactly.")


if __name__ == "__main__":
    main()
