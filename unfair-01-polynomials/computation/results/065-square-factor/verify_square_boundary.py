"""Check the seven-root families, confluent moment, and triple-cluster lift."""

import argparse
from collections import Counter
import json
from pathlib import Path

import sympy as sp

Z=sp.Symbol("z")
CASES=(
    (163,124,[1,6,121,106,68,115,117],162),
    (1279,1133,[1,1111,1202,215,1193,762,576],1277),
)


def finite(expr,p):
    return sp.Poly(expr,Z,modulus=p)


def lift(poly):
    return sp.Poly.from_list([int(c) for c in poly.all_coeffs()],Z,domain=sp.ZZ)


def hasse(poly,k,p):
    return poly.diff((0,k))*pow(int(sp.factorial(k)),-1,p)


def quotient_remainder(poly,j,p):
    rest=poly.rem(lift(j))
    assert all(int(c)%p==0 for c in rest.all_coeffs())
    return sp.Poly.from_list([int(c)//p for c in rest.all_coeffs()],Z,modulus=p)


def test_lift(g,j,k,h,p,w):
    x=g.diff() if k==3 else g
    lhs=quotient_remainder(x,j,p)
    rhs=(h*hasse(finite(g,p),k,p)*w**2).rem(j)
    return lhs==rhs


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    records=[]
    for p,h,prefix,wanted in CASES:
        counts=Counter()
        exceptional=[]
        for a in range(1,p):
            j=sp.Poly.from_list([*prefix,a],Z,modulus=p)
            invz=finite(Z,p).invert(j)
            assert int(((invz**7)*finite(Z**6,p)).rem(j).nth(6))%p==-pow(a,-1,p)%p
            collision=j.gcd(j.diff())
            if collision.degree()>0:
                exceptional.append({"constant":a,"gcd":[int(x) for x in collision.all_coeffs()]})
                assert (p,a)==(1279,144)
                factors=[(finite(Z+479,p),2),(finite(Z-308,p),1),
                         (finite(Z**2+149*Z+429,p),1),(finite(Z**2+312*Z+359,p),1)]
                assert all(g.is_irreducible for g,_ in factors)
                assert sp.prod(g**e for g,e in factors)==j and collision==finite(Z+479,p)
                continue
            counts["squarefree_parameters"]+=1
            j0=lift(j)
            w=(finite(Z**6,p)*j.diff().invert(j)).rem(j)
            b=lift(finite(h*Z**12,p).rem(j))
            base=j0**2+p*b
            assert test_lift(base,j,2,h,p,w)
            assert not test_lift(base+p,j,2,h,p,w)
            assert test_lift(base*j0,j,3,h,p,w)
            assert not test_lift((base+p)*j0,j,3,h,p,w)
            counts["lift_fixtures"]+=4
            if a==1:
                _,factors=j.factor_list()
                factor=factors[0][0]
                if factor.degree()<j.degree():
                    cofactor=j.exquo(factor)
                    mixed=base*lift(factor)
                    assert test_lift(mixed,factor,3,h,p,w.rem(factor))
                    assert test_lift(mixed,cofactor,2,h,p,w.rem(cofactor))
                    counts["mixed_cluster_fixtures"]+=2
        assert counts["squarefree_parameters"]==wanted
        assert len(exceptional)==int(p==1279)
        records.append({"prime":p,"counts":dict(counts),"exceptional_parameters":exceptional})
    report={"status":"verified","symbolic_backend":sp.__version__,"cases":records,
            "scope":"Exact synthetic polynomials test arithmetic of the necessary lift."}
    with args.out.open("x") as f:
        json.dump(report,f,indent=2)
        f.write("\n")
    print(json.dumps(report))


if __name__=="__main__":
    main()
