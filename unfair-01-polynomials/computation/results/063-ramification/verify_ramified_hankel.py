"""Rebuild the two-prime Hankel obstruction and the small collision cases."""

import argparse
from itertools import product
import json
from math import factorial, lcm
from pathlib import Path

import sympy as sp


def partitions(total, minimum=1):
    if total==0:
        yield ()
    for first in range(minimum,total+1):
        for tail in partitions(total-first,first):
            yield (first,*tail)


def permutation_orders(size):
    return {lcm(*parts) for parts in partitions(size)}


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    t=sp.Symbol("t")
    h=sp.Poly(t**6-4*t**5+6*t**4-7*t**3+4*t**2-2*t+1,t)
    assert h.discriminant()==163*1279
    assert h.eval(2)==-11 and h.eval(3)==85
    assert -11*sp.Rational(27,10)**2+8*sp.Rational(9,4)+12*sp.Rational(15,4)-12 == -sp.Rational(2919,100)
    m=1-t+t**2
    n=1-t*m
    a=(1,t,0,n,1-(1-t)*n,0,0,0)
    power_sums=[sp.S.Zero]
    for k in range(1,8):
        power_sums.append(sp.expand((-1)**(k-1)*k*a[k]+sum(
            (-1)**(j-1)*a[j]*power_sums[k-j] for j in range(1,k))))
    assert power_sums[2]==t**2 and sp.expand(power_sums[3]-t**3-3*n)==0

    matrices=[]
    cases=((163,-6,[1,157,78,152,161,51,120],111),
           (1279,168,[1,168,163,456,686,1101,875],217))
    for p,root,expected,d in cases:
        assert sp.isprime(p)
        hp=sp.Poly(h,modulus=p)
        assert hp.gcd(hp.diff())==sp.Poly(t-root,t,modulus=p)
        assert h.discriminant()%p==0 and h.discriminant()%(p*p)!=0
        assert all(factorial(k)%p for k in range(1,8))
        moments=[int(sp.diff(power_sums[k],t).subs(t,root))*pow(k,-1,p)%p for k in range(1,8)]
        assert moments==expected
        hankel=sp.Matrix(4,4,lambda i,j:moments[i+j])
        assert int(hankel.det())%p==d!=0
        matrices.append({"prime":p,"t_residue":root,"moments":moments,
                         "matrix":[list(map(int,row)) for row in hankel.tolist()],"determinant_mod_p":d})

    profiles=[]
    for deficit in range(4):
        for part in partitions(deficit):
            clusters=tuple(v+1 for v in part)
            assert all(p>max(clusters,default=1) for p, *_ in cases)
            options=sorted({e for v in product(*(permutation_orders(c) for c in clusters))
                            if (e:=lcm(*v))%2==0})
            if 4 in options:
                assert clusters==(4,)
            if 6 in options:
                assert clusters==(2,3)
            assert set(options)<={2,4,6}
            profiles.append({"gcd_degree":deficit,"cluster_sizes":clusters,
                             "possible_even_inertia_orders":options})
    assert len(profiles)==7
    result={"status":"verified","matrix_checks":matrices,"small_collision_types":profiles,
            "symbolic_backend":sp.__version__}
    with args.out.open("x") as f:
        json.dump(result,f,indent=2)
        f.write("\n")
    print(json.dumps({"status":result["status"],"determinants":[d for _,_,_,d in cases],
                      "small_collision_types":len(profiles)}))


if __name__=="__main__":
    main()
