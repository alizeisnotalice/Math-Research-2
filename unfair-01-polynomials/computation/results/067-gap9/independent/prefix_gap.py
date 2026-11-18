"""Construct the prefix equations implied by zero interior at row 9 or row 10."""

from itertools import product

import sympy as sp

T=sp.Symbol("t")


def shapes():
    return map("".join,product("ABX",repeat=3))


def pair_table(horizon,shape):
    table=dict(zip((2,3,4),shape))
    if horizon==10:table[5]="A"
    for degree in (2,3,4):
        table[horizon-degree]=shape[degree-2] if shape[degree-2]!="X" else "Z"
    table[horizon-1]="Z"
    return table


def equations(horizon,shape):
    inward=pair_table(horizon,shape)
    A={0:sp.S.One,1:T}
    B={0:sp.S.One,1:1-T}
    vars_=[]
    equations=[]
    for k,side in sorted(inward.items()):
        # Match actual indices from the current prefix rather than the positional table.
        cross=sum(A[i]*B[k-i] for i in A if i>0 and k-i in B and k-i>0)
        total=sp.expand(1-cross)
        if side=="A":
            A[k],B[k]=total,sp.S.Zero
        elif side=="B":
            A[k],B[k]=sp.S.Zero,total
        elif side=="X":
            v=sp.Symbol(f"x{k}")
            vars_.append(v)
            A[k],B[k]=v,sp.expand(total-v)
        else:
            assert side=="Z"
            A[k]=B[k]=sp.S.Zero
            equations.append(total)
    return (*vars_,T),equations,A,B


def strip(poly,symbols):
    f=sp.Poly(poly,*symbols,domain=sp.QQ)
    for t in (T,T-1):
        a=sp.Poly(t,*symbols,domain=sp.QQ)
        while f and f.rem(a).is_zero:
            f=f.exquo(a)
    return f.as_expr()


def primitive(f):
    _,p=sp.Poly(f,T,domain=sp.QQ).clear_denoms(convert=True)
    p=p.primitive()[1]
    return -p if p.LC()<0 else p
