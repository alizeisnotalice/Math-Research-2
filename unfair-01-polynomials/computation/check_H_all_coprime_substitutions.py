"""Exact integer audit for the carry formula in H(x^a)H(x^b).

This checks the finite coordinate arithmetic for 2<=a<=6, 7<=b<=12,
gcd(a,b)=1.  It does not replace the real-coefficient proof in the note.
"""
from collections import defaultdict
import argparse
import json
from math import gcd

S = (0, 1, 2, 5, 6)


def target_pairs(a, b):
    return {(i, j): a * i + b * j for i in S for j in S}


def actual_extras(a, b):
    targets = target_pairs(a, b)
    extras = {}
    for I in range(13):
        for J in range(13):
            n = a * I + b * J
            represented = [(i, j) for (i, j), value in targets.items() if value == n]
            if represented and (I, J) not in targets:
                extras[I, J] = tuple(represented)
    return extras


def expected_extras(a, b):
    out = {}
    carry_rows = {}
    allowed_i = tuple(i for i in S if i <= 12 - b)
    for j in S:
        if j >= a:
            coords = []
            for i in allowed_i:
                coord = (i + b, j - a)
                out[coord] = ((i, j),)
                coords.append(coord)
            carry_rows[j] = tuple(coords)
    return out, carry_rows


def exact_row_audit(a, b):
    extras = actual_extras(a, b)
    expected, carry_rows = expected_extras(a, b)
    assert extras == expected, (a, b, extras, expected)

    # Build the formal rows from (3).  A term is (sign, kind, exponent),
    # where H means a base H(u) contribution and D_j represents D_j(u).
    rows = defaultdict(list)
    for j in S:
        rows[j].append((+1, "H", None))
    for j in S:
        if j >= a:
            rows[j - a].append((+1, f"u^{b}D_{j}", None))
            rows[j].append((-1, f"D_{j}", None))

    # Direct coordinate computation of each formal carry: its positive and
    # negative coordinates have equal aI+bJ exactly, and the positive u
    # exponent stays <=12 exactly for the listed D_j coefficients.
    allowed_i = tuple(i for i in S if i <= 12 - b)
    for j in S:
        if j < a:
            continue
        for i in allowed_i:
            assert a * (i + b) + b * (j - a) == a * i + b * j
            assert 0 <= i + b <= 12 and 0 <= j - a <= 12

    # The degree-six row retains u^6: no D_j contains u^6 since b>=7.
    assert 6 not in allowed_i
    assert any(kind == "H" for sign, kind, _ in rows[6] if sign == 1)

    # Canonical factor row degrees lie in S.  Positive degrees summing to 6
    # have exactly the two ordered possibilities used in the proof.
    positive_degree_pairs = [(r, s) for r in S for s in S if r > 0 and s > 0 and r + s == 6]
    assert positive_degree_pairs == [(1, 5), (5, 1)]

    return {
        "a": a,
        "b": b,
        "allowed_D_u_exponents": list(allowed_i),
        "carry_rows": {str(j): [list(c) for c in coords] for j, coords in carry_rows.items()},
        "formal_v_rows": {
            str(r): [f"{sign:+d}{kind}" for sign, kind, _ in rows[r]] for r in range(7)
        },
    }


def main(output_path=None):
    rows = []
    for a in range(2, 7):
        for b in range(7, 13):
            if gcd(a, b) == 1:
                rows.append(exact_row_audit(a, b))

    # Analytic b>=13 boundary: a noncanonical q=1 coordinate would have
    # I=i+b>=b>12; q>=2 is still larger and q<0 makes I<0.
    for a in range(2, 7):
        for b in (13, 17, 19):
            if gcd(a, b) == 1:
                assert actual_extras(a, b) == {}

    rendered = json.dumps({
        "status": "passed_exact_integer_coordinate_and_row_audit",
        "S": list(S),
        "tested_pairs": len(rows),
        "rows": rows,
        "b_ge_13_argument": "q=1 gives I=i+b>12; q>=2 is larger; q<0 gives I<0",
    }, ensure_ascii=False, sort_keys=True)
    print(rendered)
    if output_path:
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write(rendered + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="optional JSON output path")
    args = parser.parse_args()
    main(args.output)
