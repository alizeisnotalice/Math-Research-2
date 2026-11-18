"""Exact tests for the dual-end seven-residue polynomials and first lift."""

import argparse
from collections import Counter
import json
from pathlib import Path

import sympy as sp


Z=sp.Symbol("z")


def field(poly,p):
    return sp.Poly(poly,Z,modulus=p)


def zz(poly):
    return sp.Poly.from_list([int(x) for x in poly.all_coeffs()],Z,domain=sp.ZZ)


def divided_remainder(g,j,p):
    remainder=g.rem(zz(j))
    assert all(int(x)%p==0 for x in remainder.all_coeffs())
    return sp.Poly.from_list([int(x)//p for x in remainder.all_coeffs()],Z,modulus=p)


def fixture(g,j,h,p):
    slope=(field(g,p).diff().diff()*pow(2,-1,p)).rem(j)
    weight=(field(Z**6,p)*j.diff().invert(j)).rem(j)
    return (h*slope*weight**2).rem(j)


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    t=sp.Symbol("t")
    h=sp.Poly(t**6-4*t**5+6*t**4-7*t**3+4*t**2-2*t+1,t)
    assert h.discriminant()==163*1279
    n=1-t+t*t-t**3
    a=[1,t,0,n,1-(1-t)*n,0,0,0]
    sums=[sp.S.Zero]
    for k in range(1,8):
        sums.append(sp.expand((-1)**(k-1)*k*a[k]+sum(
            (-1)**(j-1)*a[j]*sums[k-j] for j in range(1,k))))
    data=[]
    for p,root,lift,expected,count in (
        (163,-6,124,[1,6,121,106,68,115,117],162),
        (1279,168,1133,[1,1111,1202,215,1193,762,576],1277),
    ):
        assert h.eval(root)%p==h.diff().eval(root)%p==0
        assert (int(h.eval(root))//p)*pow(int(h.diff().diff().eval(root))//2,-1,p)%p==lift
        positive=[int(sp.diff(sums[k],t).subs(t,root))*pow(k,-1,p)%p for k in range(1,8)]
        moments=dict(zip(range(-6,7),[0]*6+positive,strict=True))
        needed=[1]
        for i in range(1,7):
            needed.append(-sum(moments[j]*needed[i-j] for j in range(1,i+1))%p)
        assert needed==expected
        matrix=sp.Matrix(7,7,lambda i,j:moments[i+j-6])
        assert int(matrix.det())%p==p-1
        totals=Counter()
        for c0 in range(p):
            j=sp.Poly.from_list([*needed,c0],Z,modulus=p)
            assert all((sum(int(j.nth(l))*moments[k+l] for l in range(8)))%p==0
                       for k in range(-6,0))
            if c0==0 or j.gcd(j.diff()).degree()>0:
                totals["ineligible"]+=1
                continue
            assert field(Z**6,p).gcd(j).degree()==0
            good=zz(j)**2+p*zz(field(lift*Z**12,p).rem(j))
            assert divided_remainder(good,j,p)==fixture(good,j,lift,p)
            bad=good+p
            assert divided_remainder(bad,j,p)!=fixture(bad,j,lift,p)
            totals["eligible"]+=1
            totals["lift_fixtures"]+=2
        assert totals["eligible"]==count
        data.append({"prime":p,"H_lift":lift,"moments_0_through_6":positive,
                     "monic_prefix":needed,"signed_hankel_determinant_mod_p":p-1,
                     "counts":dict(totals)})
    report={"status":"verified","arith_backend":sp.__version__,"parameters":data,
            "fixture_scope":"Synthetic integer polynomials for the quotient-ring lift, not binary targets."}
    with args.out.open("x") as f:
        json.dump(report,f,indent=2)
        f.write("\n")
    print(json.dumps(report))


if __name__=="__main__":
    main()
