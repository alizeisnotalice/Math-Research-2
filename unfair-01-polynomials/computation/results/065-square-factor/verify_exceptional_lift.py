"""Recompute the lift on the sole repeated family parameter at 1279."""
from pathlib import Path
import argparse,json
import sympy as sp

from verify_square_boundary import Z,finite,lift,hasse,quotient_remainder


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    p,h,r=1279,1133,-479
    j=sp.Poly.from_list([1,1111,1202,215,1193,762,576,144],Z,modulus=p)
    ell=finite(Z-r,p)
    assert j.gcd(j.diff())==ell
    other=j.exquo(ell**2)
    assert other.gcd(ell).degree()==other.gcd(other.diff()).degree()==0
    b=lift(finite(h*Z**12,p).rem(j))
    base=lift(j)**2+p*b
    rows=[]
    for c in (4,5):
        g=base*lift(ell)**(c-4)
        m=hasse(j,2,p).rem(ell)
        w=(finite(Z**6,p)*m.invert(ell)).rem(ell)
        qp=quotient_remainder(g.diff((0,c-4)),ell,p)
        check=(h*hasse(finite(g,p),c,p)*w**2).rem(ell)
        assert qp==check and not qp.is_zero
        changed=g+p*lift(ell)**(c-4)
        assert quotient_remainder(changed.diff((0,c-4)),ell,p)!=check
        ordinary=other
        w_other=(finite(Z**6,p)*j.diff().invert(ordinary)).rem(ordinary)
        assert quotient_remainder(g,ordinary,p)==(h*hasse(finite(g,p),2,p)*w_other**2).rem(ordinary)
        rows.append({"cluster_size":c,"exact_residue":int(qp.nth(0)),
                     "negative_fixture_rejected":True,"other_five_roots_checked":True})
    with args.out.open("x") as output:
        json.dump({"status":"verified","prime":p,"constant_parameter":144,"rows":rows},output,indent=2)
        output.write("\n")



if __name__=="__main__":
    main()
