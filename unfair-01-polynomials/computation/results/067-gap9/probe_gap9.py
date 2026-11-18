"""Exact prefix elimination for a first gap at 9 after fractional linear terms."""

import argparse
from itertools import product
import json
from pathlib import Path
import time

import sympy as sp


def strip_endpoint_factors(f, symbols, t):
    p = sp.Poly(f, *symbols, domain=sp.QQ)
    for endpoint in (t, t-1):
        e = sp.Poly(endpoint, *symbols, domain=sp.QQ)
        while p and p.rem(e).is_zero:
            p = p.exquo(e)
    return p.as_expr()


def model(shape):
    t = sp.Symbol("t")
    used = []
    A, B = [1,t], [1,1-t]
    zero_rows = []
    for k in range(2,9):
        value = sp.expand(1-sum(A[i]*B[k-i] for i in range(1,k)))
        if k<=4:
            row=shape[k-2]
            if row=="X":
                u=sp.Symbol(f"a{k}")
                used.append(u)
                A.append(u)
                B.append(sp.expand(value-u))
                continue
        elif k==8:
            row="X"
        else:
            row=shape[7-k]
        if row=="A":
            A.append(value); B.append(0)
        elif row=="B":
            A.append(0); B.append(value)
        else:
            A.append(0); B.append(0)
            zero_rows.append(value)
    symbols=(*used,t)
    return symbols,[strip_endpoint_factors(f,symbols,t) for f in zero_rows]


def probe(shape):
    symbols,fs=model(shape)
    start=time.monotonic()
    G=sp.groebner(fs,*symbols,domain=sp.QQ,order="lex")
    t=symbols[-1]
    eliminants=[f for f in G if f.free_symbols <= {t}]
    assert len(eliminants)<=1
    factors=[]
    for f in eliminants:
        for factor,power in sp.factor_list(f)[1]:
            p=sp.Poly(factor,t,domain=sp.ZZ).primitive()[1]
            roots=p.count_roots(0,1)
            factors.append(dict(coefficients=list(map(int,p.all_coeffs())),
                                degree=p.degree(),multiplicity=power,roots_closed_0_1=int(roots)))
    return dict(shape=shape,variables=list(map(str,symbols)),zero_rows=list(map(str,fs)),
                basis_degrees=[sp.Poly(f,*symbols).total_degree() for f in G],
                unit=(list(G)==[1]),eliminant=bool(eliminants),factors=factors,
                seconds=time.monotonic()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    with args.out.open("x") as stream:
        for choices in product("ABX",repeat=3):
            shape="".join(choices)
            if shape[:2] == "XX" or shape[0] == shape[2] == "X":
                continue
            row=probe(shape)
            stream.write(json.dumps(row)+"\n")
            stream.flush()
            print(shape,row["unit"],[(f["degree"],f["coefficients"][0],f["roots_closed_0_1"]) for f in row["factors"]],flush=True)


if __name__=="__main__":
    main()
