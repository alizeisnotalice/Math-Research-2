"""Check degree-six field invariants and the target ramification exclusions."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    x=sp.Symbol("x")
    h=sp.Poly(x**6-4*x**5+6*x**4-7*x**3+4*x**2-2*x+1,x)
    assert h.discriminant()==208477==163*1279
    assert sp.isprime(163) and sp.isprime(1279)
    signatures={}
    for p, expected in ((2,[6]),(11,[1,5]),(37,[1,2,3])):
        residue=sp.Poly(h,modulus=p)
        _, factors=residue.factor_list()
        assert all(e==1 and f.is_irreducible for f,e in factors)
        assert sp.prod(f for f,_ in factors)==residue
        signature=sorted(f.degree() for f,_ in factors)
        assert signature==expected
        signatures[p]=signature
    raw=args.packet.read_bytes()
    doc=json.loads(raw)
    rows=[]
    for task in doc["tasks"]:
        c=sp.Poly.from_list(task["C"][::-1],x,domain=sp.ZZ)
        assert c.is_monic and c.degree()==56
        radical=c.exquo(c.gcd(c.diff()))
        assert radical.is_monic and radical.gcd(radical.diff()).degree()==0
        values={}
        for p in (163,1279):
            f=sp.Poly(radical,modulus=p)
            g=f.gcd(f.diff())
            # Extended gcd independently checks the unit ideal when degree is zero.
            if g.degree()==0:
                u,v,d=f.gcdex(f.diff())
                assert d==sp.Poly(1,x,modulus=p) and u*f+v*f.diff()==d
            values[p]=g.degree()
        rows.append({"task":task["id"],"gcd_degrees":values,
                     "excluded":not all(values.values())})
    report={"status":"verified","H_discriminant":208477,"H_galois_cycle_degrees":signatures,
            "packet_sha256":sha256(raw).hexdigest(),"rows":rows,
            "excluded":sum(row["excluded"] for row in rows)}
    with args.out.open("x") as out:
        json.dump(report,out,indent=2)
        out.write("\n")
    print(json.dumps({"status":report["status"],"excluded":report["excluded"]}))


if __name__=="__main__":
    main()
