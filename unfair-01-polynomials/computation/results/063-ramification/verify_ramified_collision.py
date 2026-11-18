"""Verify the two-root local obstruction on squarefree integer targets."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


def eval_integer(coefficients, point):
    value=0
    for c in coefficients:
        value=value*point+c
    return value


def scaled_value_and_curvature(poly, root, prime):
    h=poly.all_coeffs()
    n=poly.degree()
    q=eval_integer(h,root)
    assert q%prime==0
    curvature=eval_integer([a*(n-j)*(n-j-1)//2 for j,a in enumerate(h[:-2])],root)
    assert curvature%prime!=0
    return (int(q)//prime)*pow(int(curvature),-1,prime)%prime


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    opts=parser.parse_args()
    x=sp.Symbol("x")
    h=sp.Poly(x**6-4*x**5+6*x**4-7*x**3+4*x**2-2*x+1,x,domain=sp.ZZ)
    assert h.discriminant()==163*1279
    cases=((163,-6,-2,124),(1279,168,2,1133))
    raw=opts.packet.read_bytes()
    inputs=json.loads(raw)
    rows=[]
    for task in inputs["tasks"]:
        c=sp.Poly.from_list(task["C"][::-1],x,domain=sp.ZZ)
        assert c.is_monic and c.nth(0)==1
        assert c.gcd(c.diff()).degree()==0
        invariants={}
        for p,hroot,croot,expect in cases:
            assert sp.isprime(p) and croot%p
            h_p,c_p=sp.Poly(h,modulus=p),sp.Poly(c,modulus=p)
            assert h_p.gcd(h_p.diff())==sp.Poly(x-hroot,x,modulus=p)
            assert c_p.gcd(c_p.diff())==sp.Poly(x-croot,x,modulus=p)
            left=scaled_value_and_curvature(h,hroot,p)
            assert left==expect
            right=scaled_value_and_curvature(c,croot,p)*pow(croot,-4,p)%p
            invariants[p]={"product_root":croot,"h_invariant":left,"product_invariant":right}
        assert all(row["h_invariant"]!=row["product_invariant"] for row in invariants.values())
        rows.append({"id":task["id"],"mask":task["mask"],"invariants":invariants})
    assert [row["id"] for row in rows]==list(range(12))
    result={"status":"verified","source_inputs_sha256":sha256(raw).hexdigest(),"rows":rows}
    with opts.out.open("x") as out:
        json.dump(result,out,indent=2)
        out.write("\n")
    print(json.dumps({"status":result["status"],"excluded_targets":len(rows)}))


if __name__=="__main__":
    main()
