"""Exact support-level audit for the non-binary macro-lift theorem."""


def macro_terms(a):
    return [(0, 0), (0, a), (0, 2*a), (a, -a), (a, 0)]


def audit(a, b, d):
    c = {}
    for r, e in macro_terms(a):
        for s, f in macro_terms(b):
            c[r+s, e+f] = c.get((r+s, e+f), 0) + 1
    target = set(c)
    phi = {r*d + e: (r, e) for r, e in target}
    assert len(phi) == len(target)
    actual = {(r+s, e+f) for r, e in target for s, f in target}
    for r, e in actual:
        n = r*d + e
        if n in phi:
            assert phi[n] == (r, e), (a, b, d, n, phi[n], (r, e))
    return {"a": a, "b": b, "d": d, "macro_terms": len(target),
            "max_macro_coefficient": max(c.values())}


def main():
    rows = []
    for a in range(1, 7):
        for b in range(1, 7):
            rows.append(audit(a, b, 5*(a+b) + 1))
    assert audit(1, 2, 16)["max_macro_coefficient"] > 1
    assert audit(3, 3, 31)["max_macro_coefficient"] > 1
    print({"status": "passed_exact_macro_support_v2_audit", "tested": len(rows),
           "included_pairs": ["(1,2)", "(3,3)"]})


if __name__ == "__main__":
    main()
