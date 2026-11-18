#!/usr/bin/env python3
"""Independent exact local constants for the square-factor boundary."""

import argparse
import json
from pathlib import Path

import sympy as sp


def check_lifts():
    y,p,z=sp.symbols("y p z")
    v=sp.symbols("v0:6")
    unit=sum(v[i]*y**i for i in range(6))
    for mu in (1,2):
        b=sp.symbols(f"b0:{2*mu}")
        Q=y**(2*mu)+p*sum(b[i]*y**i for i in range(2*mu))
        for fixed in (0,1):
            G=sp.Poly(Q*unit*(y+p*z)**fixed,y)
            c=2*mu+fixed
            assert sp.expand(G.nth(fixed)-p*b[0]*v[0]).coeff(p,0)==0
            assert sp.expand(G.nth(fixed)-p*b[0]*v[0]).coeff(p,1)==0
            assert G.nth(c).subs(p,0)==v[0]
    return True


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    assert not args.out.exists()
    x=sp.symbols("x")
    H=sp.Poly(x**6-4*x**5+6*x**4-7*x**3+4*x**2-2*x+1,x)
    assert sp.discriminant(H)==163*1279
    a=[sp.Integer(1),x,sp.Integer(0),1-x+x*x-x**3,1-(1-x)*(1-x+x*x-x**3)]+[sp.Integer(0)]*3
    sums=[sp.Integer(0)]
    for k in range(1,8):
        sums.append(sp.expand(sum((-1)**(i-1)*a[i]*sums[k-i] for i in range(1,k))
                              +(-1)**(k-1)*k*a[k]))
    info={}
    for p,rho in ((163,-6),(1279,168)):
        value=int(H.eval(rho))
        h2=int(H.diff().diff().eval(rho))//2
        assert value%p==0 and value%(p*p)!=0 and h2%p
        hp=value//p*pow(h2,-1,p)%p
        assert hp=={163:124,1279:1133}[p]
        moments=[int(s.diff(x).subs(x,rho))*pow(k,-1,p)%p for k,s in enumerate(sums) if k]
        assert moments[0]==1
        coefficients=[1]
        for k in range(1,7):
            coefficients.append(-sum(coefficients[k-i]*moments[i] for i in range(1,k+1))%p)
        repeats=[]
        for c in range(1,p):
            J=sp.Poly.from_list([*coefficients,c],x,modulus=p)
            g=sp.gcd(J,J.diff())
            if g.degree():
                repeats.append(c)
            # The inverse representative has leading coefficient -1/J(0).
            inverse=sp.Poly.from_list(coefficients,x,modulus=p).mul_ground(-pow(c,-1,p))
            assert (inverse*sp.Poly(x,x,modulus=p)).rem(J).is_one
            assert int(inverse.nth(6))%p==-pow(c,-1,p)%p
        assert repeats==([] if p==163 else [144])
        row=dict(lower_collision_rho=rho,h=hp,repeated_J_constants=repeats,
                 fixed_J_coefficients=coefficients)
        if p==1279:
            r=-479
            J=sp.Poly.from_list([*coefficients,144],x,modulus=p)
            Q=sp.Poly((x-r)**2,x,modulus=p)
            K,rem=divmod(J,Q)
            assert rem.is_zero and K.eval(r)!=0
            factors=[sp.Poly(f,x,modulus=p) for f in (x-308,x*x+149*x+429,x*x+312*x+359)]
            assert factors[0]*factors[1]*factors[2]==K
            assert all(f.is_irreducible for f in factors)
            beta=pow(r,6,p)*pow(int(K.eval(r)),-1,p)%p
            row.update(special_residue=r%p,K_at_r=int(K.eval(r))%p,beta=beta,
                       h_beta_square=hp*beta*beta%p,minus_inverse_J_constant=-pow(144,-1,p)%p)
        info[p]=row
    output=dict(status="verified",sympy=sp.__version__,primes=info,
                lift_derivative_identity=check_lifts())
    args.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    print(json.dumps(output,sort_keys=True))


if __name__=="__main__":
    main()
