"""Prove the impossibility of a ninth-coefficient gap by exact prefix elimination."""

import argparse
from functools import reduce
from itertools import combinations, product
import json
from pathlib import Path

import sympy as sp


def strip_endpoints(poly, symbols):
    """Cancel powers of t and t-1, which are nonzero at admissible roots."""
    t = symbols[-1]
    original = sp.Poly(poly, *symbols, domain=sp.QQ)
    for endpoint in (t, t-1):
        factor = sp.Poly(endpoint, *symbols, domain=sp.QQ)
        while original and original.rem(factor).is_zero:
            original = original.exquo(factor)
    return original.as_expr()


def equations(shape):
    t = sp.Symbol("t")
    a, b = [1, t], [1, 1-t]
    variables, fs = [], []
    for k in range(2, 9):
        budget = sp.expand(1-sum(a[i]*b[k-i] for i in range(1, k)))
        side = shape[k-2] if k <= 4 else "-" if k == 8 else shape[7-k]
        if side == "X" and k <= 4:
            value = sp.Symbol(f"a{k}")
            variables.append(value)
            a.append(value)
            b.append(sp.expand(budget-value))
        elif side == "A":
            a.append(budget)
            b.append(0)
        elif side == "B":
            a.append(0)
            b.append(budget)
        else:
            a.append(0)
            b.append(0)
            fs.append(budget)
    symbols = (*variables, t)
    return symbols, [strip_endpoints(f, symbols) for f in fs]


def eliminate(polys, var):
    retained = [f for f in polys if f and var not in f.free_symbols]
    for f, g in combinations(polys, 2):
        if f and g and var in (f.free_symbols | g.free_symbols):
            r = sp.resultant(f, g, var)
            if r:
                retained.append(r)
    return retained


def interval_horner(poly, left, right):
    lo = hi = sp.S.Zero
    for c in poly.all_coeffs():
        products = [lo*left, lo*right, hi*left, hi*right]
        lo, hi = c + min(products), c + max(products)
    return lo, hi


def isolate_negative_coefficient(shape, factor):
    # Exchange of the two factors maps BXA to AXB.
    assert shape == "AXB" or shape == "BXA"
    (u, t), (f, g) = equations(shape)
    f, g = sp.Poly(f, u), sp.Poly(g, u)
    assert f.degree() == g.degree() == 2
    affine = sp.Poly(sp.expand(g.LC()*f.as_expr()-f.LC()*g.as_expr()), u)
    assert affine.degree() == 1
    slope, constant = [sp.Poly(c, t, domain=sp.QQ) for c in affine.all_coeffs()]
    field_minpoly = factor.set_domain(sp.QQ)
    assert slope.gcd(field_minpoly).is_one
    u_value = (-constant*slope.invert(field_minpoly)).rem(field_minpoly)
    roots = factor.intervals(eps=sp.Rational(1,10**15),inf=0,sup=1)
    assert len(roots) == 1 and roots[0][1] == 1
    (left,right),_=roots[0]
    assert 0 < left < right < 1
    lo,hi=interval_horner(u_value,left,right)
    b3budget=sp.Poly(1-(1-t)*(1-t+t**2),t,domain=sp.QQ)
    if shape=="AXB":
        coefficient="A3"
    else:
        coefficient="B3"
        # In BXA row 2 is on the right, so A1*B2 is its interior term.
        b3budget=sp.Poly(1-t*(1-t+t**2),t,domain=sp.QQ)
        lo,hi=interval_horner(b3budget-u_value,left,right)
    assert hi < -sp.Rational(2,5)
    return dict(coefficient=coefficient, isolated_root=list(map(str,(left,right))),
                bound=list(map(str,(lo,hi))), primitive_affine=[str(slope.as_expr()),str(constant.as_expr())],
                reconstructed_coefficients_ascending=list(map(str,u_value.all_coeffs()[::-1])))


def check(shape):
    if shape[:2] == "XX":
        return dict(shape=shape, reason="eighth_row_at_most_one_product")
    if shape[0] == shape[2] == "X":
        return dict(shape=shape, reason="eighth_row_quadratic_upper_bound")
    symbols, fs = equations(shape)
    reduced = fs
    ranks = []
    for var in symbols[:-1]:
        reduced = eliminate(reduced, var)
        ranks.append(len(reduced))
    assert reduced
    t = symbols[-1]
    eliminant = reduce(sp.gcd, reduced)
    assert eliminant and eliminant.free_symbols <= {t}
    factors = []
    for f, exponent in sp.factor_list(eliminant)[1]:
        _, p = sp.Poly(f,t,domain=sp.QQ).clear_denoms(convert=True)
        p=p.primitive()[1]
        assert p.is_irreducible
        if p.LC()<0: p=-p
        open_count = int(p.count_roots(0,1)) - int(p.eval(0)==0) - int(p.eval(1)==0)
        row=dict(degree=int(p.degree()), coefficients=list(map(int,p.all_coeffs())), exponent=int(exponent), open_roots=open_count)
        if p.LC() != 1:
            row["reason"]="not_algebraic_integer"
        elif open_count == 0:
            row["reason"]="no_open_root"
        else:
            assert p.degree()==16 and open_count==1
            row["reason"]="negative_original_coefficient"
            row["sign"]=isolate_negative_coefficient(shape,p)
        factors.append(row)
    return dict(shape=shape, reason="elimination", elimination_counts=ranks,
                assumptions=list(map(str,fs)),factors=factors)


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    rows=[check("".join(shape)) for shape in product("ABX",repeat=3)]
    assert len(rows)==27
    exceptional=[r["shape"] for r in rows if any("sign" in f for f in r.get("factors",[]))]
    assert exceptional==["AXB","BXA"]
    report=dict(status="verified", shapes=len(rows), analytic=5, eliminations=22,
                exceptional=exceptional, sympy=sp.__version__, ledger=rows)
    with args.out.open("x") as handle:
        json.dump(report,handle,indent=2)
        handle.write("\n")
    print(report["status"],len(rows),exceptional)


if __name__=="__main__":
    main()
